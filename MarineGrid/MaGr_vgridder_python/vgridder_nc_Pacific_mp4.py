from __future__ import annotations

import glob
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import label as cc_label

import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, Point

# Integration: Using Path for rapid spatial filtering
try:
    from matplotlib.path import Path
except ImportError:
    Path = None

PROGRESS_EVERY_J = 109
EPS_RAY = 0.000002  # Match Fortran 'del' value


def nint_pos(x: float) -> int:
    return int(np.floor(x + 0.5))


def _strip_comment(line: str) -> str:
    if "<" in line:
        line = line.split("<", 1)[0]
    return line.strip()


def wrap_lon_to_360(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x < 0.0, x + 360.0, x)


def wrap_lon_to_180(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x > 180.0, x - 360.0, x)


@dataclass
class InConfig:
    title: str
    alat1: float
    alat2: float
    alon1: float
    alon2: float
    delx: float
    dely: float
    ielim_ponds: int
    nadjust: int
    layers: int
    nobarr1: int
    nobarr2: int
    numnulls: int
    nblock: int


def parse_in_file(path: str) -> InConfig:
    with open(path, "rt") as f:
        title = f.readline().rstrip("\n")

        l2 = _strip_comment(f.readline()).split()
        if len(l2) < 5 or int(float(l2[0])) != 4:
            raise ValueError(f"Unsupported window line in {path}: {l2}")
        alat1, alat2, alon1, alon2 = map(float, l2[1:5])

        l3 = _strip_comment(f.readline()).split()
        delx, dely = float(l3[0]), float(l3[1])

        ielim_ponds = int(float(_strip_comment(f.readline()).split()[0]))
        nadjust = int(float(_strip_comment(f.readline()).split()[0]))

        if nadjust > 0:
            got = 0
            while got < nadjust:
                parts = _strip_comment(f.readline()).split()
                if not parts:
                    continue
                na = int(float(parts[0]))
                got += na

        layers = int(float(_strip_comment(f.readline()).split()[0]))

        nb = _strip_comment(f.readline()).split()
        nobarr1, nobarr2 = int(float(nb[0])), int(float(nb[1]))

        _ = f.readline()

        nn = _strip_comment(f.readline()).split()
        numnulls = int(float(nn[0])) if nn else 0
        nblock = int(float(nn[1])) if len(nn) > 1 else 0

        for _ in range(1 + 1 + 1 + 1 + 5):
            f.readline()

    return InConfig(
        title=title,
        alat1=alat1,
        alat2=alat2,
        alon1=alon1,
        alon2=alon2,
        delx=delx,
        dely=dely,
        ielim_ponds=ielim_ponds,
        nadjust=nadjust,
        layers=layers,
        nobarr1=nobarr1,
        nobarr2=nobarr2,
        numnulls=numnulls,
        nblock=nblock,
    )


def build_grid(cfg: InConfig) -> Tuple[np.ndarray, np.ndarray, int, int, float, float]:
    alon0 = min(cfg.alon1, cfg.alon2)
    alon3 = max(cfg.alon1, cfg.alon2)
    alat0 = min(cfg.alat1, cfg.alat2)
    alat3 = max(cfg.alat1, cfg.alat2)

    imax = 1 + nint_pos((alon3 - alon0) / cfg.delx)
    jmax = 1 + nint_pos((alat3 - alat0) / cfg.dely)

    alon3 = alon0 + cfg.delx * (imax - 1)
    alat3 = alat0 + cfg.dely * (jmax - 1)
    if alon3 < max(cfg.alon1, cfg.alon2):
        imax += 1
    if alat3 < max(cfg.alat1, cfg.alat2):
        jmax += 1

    clon = alon0 + cfg.delx * np.arange(imax, dtype=np.float64)
    clat = alat0 + cfg.dely * np.arange(jmax, dtype=np.float64)
    return clon, clat, imax, jmax, alat0, alon0


def inside_pert_crossings(x: np.ndarray, y: np.ndarray, yi: float) -> np.ndarray:
    ys = y.astype(np.float64, copy=True)
    ys[ys == yi] += EPS_RAY

    xcr: List[float] = []
    for n in range(len(x) - 1):
        y1, y2 = ys[n], ys[n + 1]
        if yi > min(y1, y2) and yi < max(y1, y2):
            x1, x2 = float(x[n]), float(x[n + 1])
            if x1 == x2:
                xc = x1
            else:
                xc = x1 + (x2 - x1) * (yi - float(ys[n])) / (float(ys[n + 1]) - float(ys[n]))
            xcr.append(float(xc))
    return np.asarray(xcr, dtype=np.float64)


def _detect_lon_lat(col0: np.ndarray, col1: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    a = col0.astype(np.float64)
    b = col1.astype(np.float64)
    a_lat_like = np.mean(np.abs(a) <= 90.0)
    b_lat_like = np.mean(np.abs(b) <= 90.0)
    if a_lat_like > 0.95 and b_lat_like < 0.95:
        return b, a
    if b_lat_like > 0.95 and a_lat_like < 0.95:
        return a, b
    return a, b


def split_segments_by_pen(pen: np.ndarray) -> List[Tuple[int, int]]:
    bounds: List[Tuple[int, int]] = []
    isum = 0
    start = 0
    for i in range(len(pen)):
        if isum == 0:
            start = i
        isum += int(pen[i])
        if isum == 2:
            bounds.append((start, i))
            isum = 0
    if not bounds:
        bounds.append((0, len(pen) - 1))
    return bounds


def read_bp_segments(bp_path: str) -> List[Tuple[np.ndarray, np.ndarray]]:
    data = np.genfromtxt(bp_path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 3:
        raise ValueError(f"BP file has <3 cols: {bp_path}")

    col0, col1 = data[:, 0], data[:, 1]
    pen = data[:, 2].astype(np.int32)
    poly = data[:, 3].astype(np.int32) if data.shape[1] >= 4 else np.zeros(len(pen), dtype=np.int32)

    lon, lat = _detect_lon_lat(col0, col1)

    u = np.unique(pen)
    if not np.all(np.isin(u, [0, 1])):
        raise ValueError(f"[HARD FAIL] BP pen not 0/1. uniq={u[:20]} file={bp_path}")

    segs_all: List[Tuple[np.ndarray, np.ndarray, int]] = []
    for s, e in split_segments_by_pen(pen):
        pid = int(poly[s])
        segs_all.append((lon[s: e + 1].copy(), lat[s: e + 1].copy(), pid))

    main_idx = 0
    for k, (_x, _y, pid) in enumerate(segs_all):
        if pid == 1:
            main_idx = k
            break

    main_seg = (segs_all[main_idx][0], segs_all[main_idx][1])
    mask_segs = [(x, y) for i, (x, y, _pid) in enumerate(segs_all) if i != main_idx]
    segs = [main_seg] + mask_segs

    pids = sorted(set(pid for _x, _y, pid in segs_all))
    print(f"[SANITY][BP] rows={len(pen)} cols={data.shape[1]} pen_unique={u} poly_ids={pids} segments={len(segs)}")
    return segs


def read_coast_segments_by_segment_id(coast_path: str) -> List[Tuple[np.ndarray, np.ndarray]]:
    data = np.genfromtxt(coast_path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise ValueError(f"COAST file must have 4 cols (lat lon pen seg_id): {coast_path}")

    lat = data[:, 0].astype(np.float64)
    lon = data[:, 1].astype(np.float64)
    pen = data[:, 2].astype(np.int32)
    seg_id = data[:, 3].astype(np.int64)

    u = np.unique(pen)
    if not np.all(np.isin(u, [0, 1])):
        raise ValueError(f"[HARD FAIL] COAST pen not 0/1. uniq={u[:20]} file={coast_path}")

    groups: Dict[int, List[int]] = {}
    for idx, sid in enumerate(seg_id.tolist()):
        groups.setdefault(int(sid), []).append(idx)

    segs: List[Tuple[np.ndarray, np.ndarray]] = []
    for _sid, idxs in groups.items():
        xs = lon[idxs].copy()
        ys = lat[idxs].copy()
        if xs.size >= 2:
            segs.append((xs, ys))

    lengths = np.array([len(s[0]) for s in segs], dtype=np.int64)
    k = int(np.argmax(lengths))
    if k != 0:
        segs = [segs[k]] + segs[:k] + segs[k + 1:]

    print(
        f"[SANITY][COAST] rows={len(pen)} cols={data.shape[1]} pen_unique={u} "
        f"segments={len(segs)} unique_seg_id={len(groups)} longest_len={int(lengths.max())}"
    )
    return segs


def read_coast_segments_from_shapefile(shp_path: str, use_360: bool = False) -> List[Tuple[np.ndarray, np.ndarray]]:
    print(f"[READING COASTLINE] Loading shapefile: {shp_path}")
    print(f"[READING COASTLINE] Longitude mode: {'0..360' if use_360 else '-180..180'}")

    gdf = gpd.read_file(shp_path)
    segs: List[Tuple[np.ndarray, np.ndarray]] = []

    def fix_lon(lon: np.ndarray) -> np.ndarray:
        lon = np.asarray(lon, dtype=np.float64)
        return wrap_lon_to_360(lon) if use_360 else wrap_lon_to_180(lon)

    def extract_coords(geom):
        if isinstance(geom, LineString):
            coords = np.array(geom.coords, dtype=np.float64)
            lon = fix_lon(coords[:, 0])
            lat = coords[:, 1]
            segs.append((lon, lat))
        elif isinstance(geom, MultiLineString):
            for line in geom.geoms:
                coords = np.array(line.coords, dtype=np.float64)
                lon = fix_lon(coords[:, 0])
                lat = coords[:, 1]
                segs.append((lon, lat))
        elif isinstance(geom, Polygon):
            coords = np.array(geom.exterior.coords, dtype=np.float64)
            lon = fix_lon(coords[:, 0])
            lat = coords[:, 1]
            segs.append((lon, lat))
        elif isinstance(geom, MultiPolygon):
            for poly in geom.geoms:
                coords = np.array(poly.exterior.coords, dtype=np.float64)
                lon = fix_lon(coords[:, 0])
                lat = coords[:, 1]
                segs.append((lon, lat))

    for geom in gdf.geometry:
        if geom is not None:
            extract_coords(geom)

    if segs:
        lengths = np.array([len(s[0]) for s in segs], dtype=np.int64)
        k = int(np.argmax(lengths))
        if k != 0:
            segs = [segs[k]] + segs[:k] + segs[k + 1:]

    print(f"[SANITY][COAST_SHP] Extracted {len(segs)} segments. Longest segment length={len(segs[0][0]) if segs else 0}")
    if segs:
        xs_all = np.concatenate([s[0] for s in segs])
        print(f"[SANITY][COAST_SHP] lon range after wrap: {xs_all.min()} to {xs_all.max()}")

    return segs


def read_adcirc_full_mesh(adcirc_path: str):
    with open(adcirc_path, "rt") as f:
        _title = f.readline().rstrip("\n")
        parts = f.readline().split()
        nt, nodes = int(float(parts[0])), int(float(parts[1]))
        xs, ys, ds = np.zeros(nodes), np.zeros(nodes), np.zeros(nodes)
        for _ in range(nodes):
            row = f.readline().split()
            nid = int(float(row[0])) - 1
            xs[nid], ys[nid], ds[nid] = float(row[1]), float(row[2]), float(row[3])
        ncellnodes = np.zeros((nt, 3), dtype=np.int32)
        for _ in range(nt):
            row = f.readline().split()
            eid = int(float(row[0])) - 1
            ncellnodes[eid, :] = [int(float(x)) - 1 for x in row[2:5]]
    return xs, ys, ds, ncellnodes


def read_adcirc_netcdf_dual_path(mesh_path: str, datum_path: str):
    print(f"[READING ADCIRC MESH] Loading: {mesh_path}")
    with Dataset(mesh_path, 'r') as ncm:
        if 'x' in ncm.variables:
            xs = ncm.variables['x'][:]
        elif 'lon' in ncm.variables:
            xs = ncm.variables['lon'][:]
        else:
            raise KeyError(f"Mesh file {mesh_path} missing longitude (checked x, lon)")

        if 'y' in ncm.variables:
            ys = ncm.variables['y'][:]
        elif 'lat' in ncm.variables:
            ys = ncm.variables['lat'][:]
        else:
            raise KeyError(f"Mesh file {mesh_path} missing latitude (checked y, lat)")

        elements = ncm.variables['element'][:] if 'element' in ncm.variables else ncm.variables['ele'][:]
        ncellnodes = np.array(elements, dtype=np.int32) - 1

    print(f"[READING ADCIRC DATUM] Loading: {datum_path}")
    with Dataset(datum_path, 'r') as ncd:
        datum_var = None
        for v in ['d', 'mllw', 'z', 'depth', 'zeta', 'datum']:
            if v in ncd.variables:
                datum_var = v
                break

        if datum_var:
            ds = ncd.variables[datum_var][:]
            print(f"  Found datum variable: {datum_var}")
        else:
            raise ValueError(f"Could not find datum variable in {datum_path}. Checked for: d, mllw, etc.")

        if np.ma.isMaskedArray(ds):
            ds = ds.filled(9.999)

    return np.array(xs).flatten(), np.array(ys).flatten(), np.array(ds).flatten(), ncellnodes


def precompute_bboxes(segs: Sequence[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    bb = np.empty((len(segs), 4), dtype=np.float64)
    for i, (x, y) in enumerate(segs):
        bb[i, 0] = float(np.min(x))
        bb[i, 1] = float(np.max(x))
        bb[i, 2] = float(np.min(y))
        bb[i, 3] = float(np.max(y))
    return bb


def prune_segments_to_grid_window(segs, bbox, lon0, lon3, lat0, lat3):
    keep = ~(
        (bbox[:, 1] < lon0 - EPS_RAY) |
        (bbox[:, 0] > lon3 + EPS_RAY) |
        (bbox[:, 3] < lat0 - EPS_RAY) |
        (bbox[:, 2] > lat3 + EPS_RAY)
    )
    idx = np.where(keep)[0]
    segs2 = [segs[i] for i in idx.tolist()]
    bbox2 = bbox[keep]
    return segs2, bbox2


def _step3_worker_chunk(args):
    rows, clon, clat, bp_segs, coast_segs, coast_bbox = args
    imax = len(clon)
    block = np.zeros((imax, len(rows)), dtype=np.int8)

    for jj, j in enumerate(rows):
        yt = float(clat[j])
        bp_xcr = [inside_pert_crossings(s[0], s[1], yt) for s in bp_segs]
        if all(x.size == 0 for x in bp_xcr):
            continue

        active_c = np.where((coast_bbox[:, 2] - EPS_RAY <= yt) & (yt <= coast_bbox[:, 3] + EPS_RAY))[0]
        c_active = [inside_pert_crossings(coast_segs[i][0], coast_segs[i][1], yt) for i in active_c]

        is_in_bp = np.zeros(imax, dtype=bool)
        xcr0 = bp_xcr[0]
        if xcr0.size > 0:
            counts0 = np.sum(clon[:, np.newaxis] < xcr0, axis=1)
            is_in_bp = (counts0 % 2 == 1)
            for xcr in bp_xcr[1:]:
                if xcr.size > 0:
                    counts_m = np.sum(clon[:, np.newaxis] < xcr, axis=1)
                    is_in_bp &= ~(counts_m % 2 == 1)

        is_land = np.zeros(imax, dtype=bool)
        for xcr in c_active:
            if xcr.size > 0:
                counts_c = np.sum(clon[:, np.newaxis] < xcr, axis=1)
                is_land |= (counts_c % 2 == 1)

        block[:, jj] = (is_in_bp & ~is_land).astype(np.int8)

    return int(rows[0]), block


def ponds_label(jfield: np.ndarray):
    wet = (jfield > 0).astype(np.int8)
    labeled, num = cc_label(wet)
    if num <= 0:
        return labeled.astype(np.int32), 0
    sizes = np.bincount(labeled.ravel())
    main_id = int(np.argmax(sizes[1:]) + 1)
    return labeled.astype(np.int32), main_id


def ponds_field_from_ipond(ipond: np.ndarray, main_id: int) -> np.ndarray:
    out = np.full(ipond.shape, -88, dtype=np.int32)
    if main_id > 0:
        mask = (ipond > 0) & (ipond != main_id)
        out[mask] = ipond[mask]
    return out


def apply_pond_removal(jfield, ipond, main_id, jset, ielim_ponds):
    if ielim_ponds <= 0 or main_id <= 0:
        return
    mask = (ipond > 0) & (ipond != main_id)
    jfield[mask] = jset


def make_layers(jfield: np.ndarray, layers: int, nblock: int) -> np.ndarray:
    if layers <= 0:
        return jfield
    imax, jmax = jfield.shape
    out = jfield.copy()

    for l in range(1, layers + 1):
        jnew = 1 + l
        ifield0 = np.zeros_like(out, dtype=np.int16)
        for i in range(1, imax - 1):
            for j in range(1, jmax - 1):
                if out[i, j] != 0:
                    continue
                neigh = (out[i - 1, j], out[i + 1, j], out[i, j - 1], out[i, j + 1])
                idx1 = int(max(neigh))
                idx2 = int(min(neigh))
                if idx1 > 0:
                    ifield0[i, j] = 1
                if nblock == 1 and idx2 < 0:
                    ifield0[i, j] = -1

        nflip = 0
        for i in range(1, imax - 1):
            for j in range(1, jmax - 1):
                if ifield0[i, j] == 0:
                    continue
                if ifield0[i, j] == 1:
                    out[i, j] = jnew
                    nflip += 1
                else:
                    out[i, j] = -jnew

        print(f"  layer={l} jnew={jnew} num added={nflip}")
        if nflip == 0:
            break

    out[out < -1] = 0
    return out


def ck_bp_after_layers(jfield, bp_segs, clon, clat):
    out = jfield.copy()
    imax, jmax = out.shape
    nch = 0
    
    print("  -> Running vectorized boundary & hole enforcement with Bounding Box pre-filters...")
    
    # PRE-CALCULATE BOUNDING BOXES FOR ALL HOLES TO SAVE MASSIVE TIME
    hole_bboxes = []
    for seg in bp_segs[1:]:
        hole_bboxes.append({
            'ymin': np.min(seg[1]),
            'ymax': np.max(seg[1]),
            'seg': seg
        })
        
    for j in range(jmax):
        # Only process the row if there are actual water/layer points in it
        row_mask = (out[:, j] > 0)
        if not np.any(row_mask):
            continue
            
        yt = float(clat[j])
        
        # 1. Enforce Main Outer Boundary
        xcr0 = inside_pert_crossings(bp_segs[0][0], bp_segs[0][1], yt)
        if xcr0.size == 0:
            out[row_mask, j] = 0
            nch += np.sum(row_mask)
            continue
            
        # Vectorized check for the main boundary
        xcr0_sorted = np.sort(xcr0)
        crossings_right = len(xcr0_sorted) - np.searchsorted(xcr0_sorted, clon, side='right')
        outside_main = (crossings_right % 2 == 0)
        invalid_mask = outside_main & row_mask
        
        # 2. Enforce the 49 Holes 
        for hole in hole_bboxes:
            # INSTANTLY SKIP the hole if the row is above or below the island!
            if yt < hole['ymin'] or yt > hole['ymax']:
                continue 
                
            xcr_h = inside_pert_crossings(hole['seg'][0], hole['seg'][1], yt)
            if xcr_h.size > 0:
                xcr_h_sorted = np.sort(xcr_h)
                cross_h_right = len(xcr_h_sorted) - np.searchsorted(xcr_h_sorted, clon, side='right')
                inside_hole = (cross_h_right % 2 == 1)
                
                # Add any points that fell inside a hole to the invalid mask
                invalid_mask |= (inside_hole & row_mask)
                
        # Bulk-delete all invalid points in this row at once!
        to_clear = invalid_mask & row_mask
        num_to_clear = np.sum(to_clear)
        if num_to_clear > 0:
            out[to_clear, j] = 0
            nch += num_to_clear

    print(f"  number changed (ck_bp)={nch}")
    return out

def remove_barriers2(jfield: np.ndarray, nobarr2: int) -> np.ndarray:
    imax, jmax = jfield.shape
    isv = np.minimum(1, jfield).astype(np.int32, copy=True)
    iter_no = 0
    ntotal = 0

    while True:
        iter_no += 1
        numc = 0

        for i in range(1, imax - 2):
            for j in range(1, jmax - 1):
                isum_left = isv[i, j] + isv[i, j + 1]
                isum_cent = isv[i + 1, j] + isv[i + 1, j + 1]
                isum_rght = isv[i + 2, j] + isv[i + 2, j + 1]
                isum_bott = isv[i, j] + isv[i + 2, j]
                isum_topp = isv[i, j + 1] + isv[i + 2, j + 1]
                icell1 = isum_left + isum_cent
                icell2 = isum_rght + isum_cent
                if icell1 == 0 or icell2 == 0 or isum_cent != 0:
                    continue
                numc += 1
                laynew = int(np.max(jfield[i:i + 3, j:j + 2])) + 1
                for k in (0, 1):
                    ns = isum_bott if k == 0 else isum_topp
                    if ns > 0 and jfield[i + 1, j + k] == 0:
                        jfield[i + 1, j + k] = laynew
                        isv[i + 1, j + k] = 1
                        ntotal += 1

        for i in range(1, imax - 1):
            for j in range(1, jmax - 2):
                isum_bott = isv[i, j] + isv[i + 1, j]
                isum_cent = isv[i, j + 1] + isv[i + 1, j + 1]
                isum_topp = isv[i, j + 2] + isv[i + 1, j + 2]
                isum_left = isv[i, j] + isv[i, j + 2]
                isum_rght = isv[i + 1, j] + isv[i + 1, j + 2]
                icell1 = isum_bott + isum_cent
                icell2 = isum_topp + isum_cent
                if icell1 == 0 or icell2 == 0 or isum_cent != 0:
                    continue
                numc += 1
                laynew = int(np.max(jfield[i:i + 2, j:j + 3])) + 1
                for k in (0, 1):
                    ns = isum_left if k == 0 else isum_rght
                    if ns > 0 and jfield[i + k, j + 1] == 0:
                        jfield[i + k, j + 1] = laynew
                        isv[i + k, j + 1] = 1
                        ntotal += 1

        if numc == 0:
            break

    print(f"  iter={iter_no} ntotal changed={ntotal}")
    return jfield


def write_nc_gtx_like(fname, field_imax_jmax, alat0, alon0, dely, delx, dtype):
    imax, jmax = field_imax_jmax.shape
    box_lon = alon0 + (360.0 if alon0 < 0 else 0.0)
    box = np.array([alat0, box_lon, dely, delx], dtype=np.float64)

    with Dataset(fname, "w", format="NETCDF4") as nc:
        nc.createDimension("box", 4)
        nc.createDimension("longitude", imax)
        nc.createDimension("latitude", jmax)

        v_box = nc.createVariable("box", "f8", ("box",))
        v_box[:] = box

        v_field = nc.createVariable("field", dtype, ("latitude", "longitude"), zlib=True, complevel=8)
        v_field[:, :] = field_imax_jmax.T

    print(f"Wrote: {os.path.basename(fname)}")


def sanity_values(tag: str, arr: np.ndarray) -> None:
    u = np.unique(arr)
    short = u if u.size <= 20 else np.concatenate([u[:10], u[-10:]])
    print(f"[SANITY] {tag}: shape={arr.shape} min={arr.min()} max={arr.max()} uniq={short}")


def inside_tri3_fortran(xtri, ytri, xp, yp):
    def r_fortran(z):
        return z if z >= 0 else z + (2 * np.pi)

    a12, a13 = np.arctan2(ytri[1] - ytri[0], xtri[1] - xtri[0]), np.arctan2(ytri[2] - ytri[0], xtri[2] - xtri[0])
    d1, d2 = r_fortran(a12 - a13), r_fortran(a13 - a12)
    d, a0 = (d2, a12) if d1 > np.pi else (d1, a13)
    if r_fortran(np.arctan2(yp - ytri[0], xp - xtri[0]) - a0) > d:
        return False

    a23, a21 = np.arctan2(ytri[2] - ytri[1], xtri[2] - xtri[1]), np.arctan2(ytri[0] - ytri[1], xtri[0] - xtri[1])
    d1, d2 = r_fortran(a23 - a21), r_fortran(a21 - a23)
    d, a0 = (d2, a23) if d1 > np.pi else (d1, a21)
    if r_fortran(np.arctan2(yp - ytri[1], xp - xtri[1]) - a0) > d:
        return False

    return True


def create_icells(n_nodes, ncellnodes):
    icells = [[] for _ in range(n_nodes)]
    for e in range(ncellnodes.shape[0]):
        for k in range(3):
            icells[int(ncellnodes[e, k])].append(e)
    return icells


def ck_nodes_v7(jfield, clon, clat, mesh_path, datum_path, numnulls):
    print("\n <ck_nodes>")
    if mesh_path.lower().endswith('.nc'):
        xs, ys, _, ncellnodes = read_adcirc_netcdf_dual_path(mesh_path, datum_path)
    else:
        xs, ys, _, ncellnodes = read_adcirc_full_mesh(mesh_path)

    with Dataset(datum_path, 'r') as ncd:
        datum_var = None
        for v in ['d', 'mllw', 'z', 'depth', 'zeta', 'datum']:
            if v in ncd.variables:
                datum_var = v
                break
        ds = ncd.variables[datum_var][:].flatten()
        if np.ma.isMaskedArray(ds):
            ds = ds.filled(9.999)

    icells = create_icells(len(xs), ncellnodes)

    imax, jmax = jfield.shape
    dx, dy = float(clon[1] - clon[0]), float(clat[1] - clat[0])
    x0, y0 = float(clon[0]), float(clat[0])

    lon_min, lon_max = clon[0], clon[-1]
    lat_min, lat_max = clat[0], clat[-1]

    idry = np.zeros_like(jfield, dtype=np.int8)
    ndry_count = 0

    for n in range(len(xs)):
        if ds[n] < 9.0:
            continue
        if not (lon_min - dx <= xs[n] <= lon_max + dx and lat_min - dy <= ys[n] <= lat_max + dy):
            continue

        ill = int(((xs[n] - x0) / dx) + 1e-9)
        jll = int(((ys[n] - y0) / dy) + 1e-9)

        if not (0 <= ill < imax - 1 and 0 <= jll < jmax - 1):
            continue

        ndry_count += 1
        for i in range(ill, ill + 2):
            for j in range(jll, jll + 2):
                if idry[i, j] == 1:
                    continue
                xp, yp = float(clon[i]), float(clat[j])
                for kc in icells[n]:
                    nod = ncellnodes[kc]
                    if inside_tri3_fortran(xs[nod], ys[nod], xp, yp):
                        num9 = sum(1 for node_idx in nod if ds[node_idx] > 9.0)
                        if num9 >= numnulls:
                            idry[i, j] = 1
                            break

    static_idry = idry.copy()
    for i in range(1, imax - 1):
        for j in range(1, jmax - 1):
            if static_idry[i, j] == 0:
                if (static_idry[i-1, j] + static_idry[i+1, j] +
                    static_idry[i, j-1] + static_idry[i, j+1]) == 4:
                    idry[i, j] = 1

    print(f"  ndry={ndry_count}")
    out = jfield.copy()
    out[idry == 1] = -1
    return out

# =============================================================================
# INTEGRATED ULTRA-FAST INDEX-BASED PRE-WETTING
# =============================================================================

def get_wet_mesh_mask(clon, clat, mesh_path, bp_coords, coast_segs, depth_limit=20.0):
    if Path is None:
        return np.zeros((len(clon), len(clat)), dtype=bool)

    print(f"[MESH] Index-based river mapping (Depth Filter: < {depth_limit}m)...")
    xs, ys, ds, ncellnodes = read_adcirc_full_mesh(mesh_path)
    
    # ---------------------------------------------------------
    # LONGITUDE WRAP FIX
    # ---------------------------------------------------------
    if np.nanmax(clon) > 180.0 and np.nanmin(xs) < 0.0:
        xs = np.where(xs < 0.0, xs + 360.0, xs)
    elif np.nanmax(clon) <= 180.0 and np.nanmax(xs) > 180.0:
        xs = np.where(xs > 180.0, xs - 360.0, xs)
        
    xmin, xmax = bp_coords[:, 0].min(), bp_coords[:, 0].max()
    ymin, ymax = bp_coords[:, 1].min(), bp_coords[:, 1].max()
    
    # =========================================================
    # BUILD THE TRUE SHORELINE FOR THE 2-METER BUFFER
    # =========================================================
    lines = [np.column_stack((x, y)) for x, y in coast_segs if len(x) >= 2]
    shoreline = MultiLineString(lines) if lines else None

    tri_xs, tri_ys = xs[ncellnodes], ys[ncellnodes]
    in_win = ~((tri_xs.max(axis=1) < xmin) | (tri_xs.min(axis=1) > xmax) | 
               (tri_ys.max(axis=1) < ymin) | (tri_ys.min(axis=1) > ymax))
    
    avg_depths = np.mean(ds[ncellnodes], axis=1)
    
    # THE FIX: Exclude the 9.999 dry/missing nodes! Water must be > 0.0 AND < 9.0
    wet_els = ncellnodes[in_win & (np.sum((ds[ncellnodes] > 0.0) & (ds[ncellnodes] < 9.0), axis=1) >= 2) & (avg_depths < depth_limit)]
    
    dx, dy = float(clon[1] - clon[0]), float(clat[1] - clat[0])
    x0, y0 = float(clon[0]), float(clat[0])
    imax, jmax = len(clon), len(clat)
    
    is_wet = np.zeros((imax, jmax), dtype=bool)
    is_burned = np.zeros((imax, jmax), dtype=bool)
    bp_path = Path(bp_coords)

    for e in wet_els:
        tx, ty = xs[e], ys[e]
        i_min, i_max = max(0, int((tx.min() - x0) / dx)), min(imax - 1, int((tx.max() - x0) / dx) + 1)
        j_min, j_max = max(0, int((ty.min() - y0) / dy)), min(jmax - 1, int((ty.max() - y0) / dy) + 1)
        
        if i_max <= i_min or j_max <= j_min: continue

        ii, jj = np.meshgrid(clon[i_min:i_max], clat[j_min:j_max], indexing='ij')
        mask = Path(np.column_stack((tx, ty))).contains_points(np.column_stack((ii.ravel(), jj.ravel()))).reshape(ii.shape)
        
        if np.any(mask):
            is_wet[i_min:i_max, j_min:j_max] |= mask
        else:
            # VIRTUAL BURN
            cx, cy = np.mean(tx), np.mean(ty)
            nearest_i = min(imax - 1, max(0, int(round((cx - x0) / dx))))
            nearest_j = min(jmax - 1, max(0, int(round((cy - y0) / dy))))
            
            dot_lon = x0 + nearest_i * dx
            dot_lat = y0 + nearest_j * dy
            
            # TRUE SHORELINE BUFFER CHECK
            if shoreline is not None:
                dist_meters = shoreline.distance(Point(dot_lon, dot_lat)) * 111111.0
                if dist_meters <= 2.0:  
                    is_wet[nearest_i, nearest_j] = True
                    is_burned[nearest_i, nearest_j] = True

    # FINAL CLEANUP STAGE 2/2
    ii, jj = np.where(is_wet & ~is_burned)
    if len(ii) > 0:
        valid = bp_path.contains_points(np.column_stack((clon[ii], clat[jj])))
        is_wet[ii[~valid], jj[~valid]] = False
        
    return is_wet


def main():
    print(" PROGRAM vgridder.py")
    ref_dir = os.environ.get("REF_DIR", "").strip()
    bp_root = os.environ.get("BP_ROOT", "").strip()
    coast_path = os.environ.get("COAST_PATH", "").strip()
    datum_path = os.environ.get("DATUM_FILE", "").strip()
    mesh_path = os.environ.get("MESH_PATH", "").strip()

    folders = sorted(glob.glob(os.path.join(bp_root, "*_01")))
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    if task_id < 0 or task_id >= len(folders):
        return

    target = folders[task_id]
    folder_name = os.path.basename(target)
    cfg = parse_in_file(os.path.join(ref_dir, f"vgridder_{folder_name}.in"))
    clon, clat, imax, jmax, alat0, alon0 = build_grid(cfg)

    print(f"\n STEP 1. Read in data\n Folder: {folder_name}\n Grid: imax={imax} jmax={jmax}")

    bp_segs = read_bp_segments(os.path.join(target, "cpolygon_xyij01.dat"))

    if coast_path.lower().endswith('.shp'):
        use_360 = float(np.nanmax(clon)) > 180.0
        coast_segs = read_coast_segments_from_shapefile(coast_path, use_360=use_360)
    else:
        coast_segs = read_coast_segments_by_segment_id(coast_path)

    coast_bbox = precompute_bboxes(coast_segs)
    coast_segs, coast_bbox = prune_segments_to_grid_window(coast_segs, coast_bbox, clon[0], clon[-1], clat[0], clat[-1])

    print("\n STEP 3: Create a Marine Grid (Vectorized + MP)")
    jfield = np.zeros((imax, jmax), dtype=np.int32)
    t_start = time.time()
    nworkers = max(1, int(os.environ.get("VGRID_NWORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1"))))

    if nworkers <= 1:
        start_j, block = _step3_worker_chunk((np.arange(jmax, dtype=int), clon, clat, bp_segs, coast_segs, coast_bbox))
        jfield[:, start_j:start_j + block.shape[1]] = block.astype(np.int32)
    else:
        nworkers = min(nworkers, jmax)
        chunks = np.array_split(np.arange(jmax, dtype=int), nworkers)
        job_args = [(chunk, clon, clat, bp_segs, coast_segs, coast_bbox) for chunk in chunks if len(chunk) > 0]
        print(f"[STEP3] Using {len(job_args)} worker processes")
        with ProcessPoolExecutor(max_workers=len(job_args)) as ex:
            for start_j, block in ex.map(_step3_worker_chunk, job_args):
                jfield[:, start_j:start_j + block.shape[1]] = block.astype(np.int32)
                done_j = start_j + block.shape[1]
                print(f"[PROGRESS] rows={done_j}/{jmax} elapsed={time.time()-t_start:.1f}s water={int(np.sum(jfield > 0))}")

    sanity_values("jfield initial", jfield)
    
    # -------------------------------------------------------------------------
    # INTEGRATED RIVER PRE-WETTING (STAGE 1/2)
    # -------------------------------------------------------------------------
    if mesh_path and os.path.exists(mesh_path):
        mesh_wet_mask = get_wet_mesh_mask(clon, clat, mesh_path, np.column_stack(bp_segs[0]), coast_segs)
        jfield[mesh_wet_mask] = 1
        river_points = mesh_wet_mask.copy() # Store for physical fix
        print(f"[MESH] Pre-wetting integrated {np.sum(mesh_wet_mask)} river points.")
    else:
        river_points = np.zeros_like(jfield, dtype=bool)

    # =====================================================================
    # THE COMMAND FILE ARCHITECTURE (LOCAL OVERRIDES)
    # =====================================================================
    override_file = os.path.join(target, "local_overrides.dat")
    if os.path.exists(override_file):
        print(f"\n[FIX] Reading manual commands from {override_file}...")
        lon_grid, lat_grid = np.meshgrid(clon, clat, indexing='ij')
        
        with open(override_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                # Skip empty lines or comments
                if len(parts) < 5 or line.startswith("#"): 
                    continue
                    
                cmd = parts[0].upper()
                lon_min, lon_max = float(parts[1]), float(parts[2])
                lat_min, lat_max = float(parts[3]), float(parts[4])
                
                local_mask = (
                    (lon_grid >= lon_min) & (lon_grid <= lon_max) & 
                    (lat_grid >= lat_min) & (lat_grid <= lat_max)
                )
                
                if cmd == "FORCE_WATER":
                    jfield[local_mask] = 1
                    river_points[local_mask] = True  # Protects from the Pond Eraser!
                    print(f"  -> {cmd}: Turned {int(np.sum(local_mask))} dots to WATER.")
                    
                elif cmd == "FORCE_LAND":
                    jfield[local_mask] = 0
                    river_points[local_mask] = False # Removes river protection
                    print(f"  -> {cmd}: Turned {int(np.sum(local_mask))} dots to LAND.")

    # =====================================================================
    # SHAPEFILE OVERRIDE ARCHITECTURE (force_water.shp & force_land.shp)
    # Uses a MASTER file in REF_DIR and instantly filters via bounding box!
    # =====================================================================
    force_water_shp = os.path.join(ref_dir, "force_water.shp")
    if os.path.exists(force_water_shp):
        print(f"\n[FIX] Reading manual WATER overrides from master file: {force_water_shp}...")
        try:
            water_segs = read_coast_segments_from_shapefile(force_water_shp, use_360=(float(np.nanmax(clon)) > 180.0))
            if water_segs and Path is not None:
                lon_grid, lat_grid = np.meshgrid(clon, clat, indexing='ij')
                grid_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
                total_forced = 0
                
                grid_min_lon, grid_max_lon = clon[0], clon[-1]
                grid_min_lat, grid_max_lat = clat[0], clat[-1]
                
                for wx, wy in water_segs:
                    # Prune polygons that don't belong to this island's domain!
                    if np.max(wx) < grid_min_lon or np.min(wx) > grid_max_lon or \
                       np.max(wy) < grid_min_lat or np.min(wy) > grid_max_lat:
                        continue
                        
                    poly_path = Path(np.column_stack((wx, wy)))
                    mask_flat = poly_path.contains_points(grid_points)
                    mask_2d = mask_flat.reshape(imax, jmax)
                    
                    if np.any(mask_2d):
                        jfield[mask_2d] = 1
                        if river_points is not None:
                            river_points[mask_2d] = True  # Protects from the Pond Eraser!
                        total_forced += np.sum(mask_2d)
                
                print(f"  -> SHAPEFILE: Turned {total_forced} dots to WATER in this domain.")
        except Exception as e:
            print(f"  -> [WARNING] Failed to process {force_water_shp}: {e}")

    force_land_shp = os.path.join(ref_dir, "force_land.shp")
    if os.path.exists(force_land_shp):
        print(f"\n[FIX] Reading manual LAND overrides from master file: {force_land_shp}...")
        try:
            land_segs = read_coast_segments_from_shapefile(force_land_shp, use_360=(float(np.nanmax(clon)) > 180.0))
            if land_segs and Path is not None:
                lon_grid, lat_grid = np.meshgrid(clon, clat, indexing='ij')
                grid_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
                total_forced = 0
                
                grid_min_lon, grid_max_lon = clon[0], clon[-1]
                grid_min_lat, grid_max_lat = clat[0], clat[-1]
                
                for lx, ly in land_segs:
                    # Prune polygons that don't belong to this island's domain!
                    if np.max(lx) < grid_min_lon or np.min(lx) > grid_max_lon or \
                       np.max(ly) < grid_min_lat or np.min(ly) > grid_max_lat:
                        continue
                        
                    poly_path = Path(np.column_stack((lx, ly)))
                    mask_flat = poly_path.contains_points(grid_points)
                    mask_2d = mask_flat.reshape(imax, jmax)
                    
                    if np.any(mask_2d):
                        jfield[mask_2d] = 0
                        if river_points is not None:
                            river_points[mask_2d] = False # Removes river protection
                        total_forced += np.sum(mask_2d)
                        
                print(f"  -> SHAPEFILE: Turned {total_forced} dots to LAND in this domain.")
        except Exception as e:
            print(f"  -> [WARNING] Failed to process {force_land_shp}: {e}")
    # =====================================================================
    write_nc_gtx_like(f"water_{folder_name}.nc", (jfield > 0).astype(np.int8), alat0, alon0, cfg.dely, cfg.delx, "i1")

    ip_init, mid_init = ponds_label(jfield)
    write_nc_gtx_like(f"ponds_{folder_name}.nc", ponds_field_from_ipond(ip_init, mid_init), alat0, alon0, cfg.dely, cfg.delx, "i4")
    apply_pond_removal(jfield, ip_init, mid_init, 0, cfg.ielim_ponds)

    if cfg.numnulls > 0 and datum_path and mesh_path:
        print(f"\n STEP 4: Check for dry nodes\n numnulls={cfg.numnulls}")
        jfield = ck_nodes_v7(jfield, clon, clat, mesh_path, datum_path, cfg.numnulls)
        sanity_values("jfield after ck_nodes", jfield)
        ip_dry, mid_dry = ponds_label(jfield)
        apply_pond_removal(jfield, ip_dry, mid_dry, -1, cfg.ielim_ponds)

    write_nc_gtx_like(f"dry_{folder_name}.nc", (jfield > 0).astype(np.int8), alat0, alon0, cfg.dely, cfg.delx, "i1")

    # -------------------------------------------------------------------------
    # THE TIMELINE SHIFT (STAGE 2/2)
    # Resurrect the rivers BEFORE layers run so the engine can see them!
    # -------------------------------------------------------------------------
    if river_points is not None:
        print("[FIX] Re-applying river water BEFORE layers flow upstream...")
        # Restore the water, but respect any nodes explicitly killed by ck_nodes (-1)
        mask_to_restore = river_points & (jfield != -1)
        jfield[mask_to_restore] = 1

    # -------------------------------------------------------------------------
    # THE TINY ISLAND ERASER (Now with Diagonal Connectivity)
    # -------------------------------------------------------------------------
    print("[FIX] Erasing tiny rogue islands...")
    land_mask = (jfield == 0).astype(np.int8)
    
    # Create a 3x3 matrix to tell the script that diagonal dots are connected
    diag_structure = np.ones((3, 3), dtype=int) 
    
    labeled_land, num_land_features = cc_label(land_mask, structure=diag_structure)
    if num_land_features > 0:
        land_sizes = np.bincount(labeled_land.ravel())
        for land_id in range(1, num_land_features + 1):
            if land_sizes[land_id] <= 1:
                jfield[labeled_land == land_id] = 1
                
    if cfg.layers > 0:
        print(f"\n STEP 6: Add layers if needed\n layers={cfg.layers}")
        jfield = make_layers(jfield, cfg.layers, cfg.nblock)
        
        sanity_values("jfield after layers + river fix", jfield)

    print("\n STEP 7: Remove points outside the BP (and inside holes)")
    jfield = ck_bp_after_layers(jfield, bp_segs, clon, clat)
    sanity_values("jfield after ck_bp", jfield)

    if cfg.nobarr1 >= 1:
        print(f"\n STEP 9: Eliminate barriers\n nobarr1={cfg.nobarr1}")
        jfield = remove_barriers2(jfield, cfg.nobarr2)
        sanity_values("jfield after remove_barriers2", jfield)

    marine_out = jfield.copy()
    marine_out[marine_out < 0] = 0
    write_nc_gtx_like(f"marine_{folder_name}.nc", marine_out.astype(np.int8), alat0, alon0, cfg.dely, cfg.delx, "i1")


if __name__ == "__main__":
    main()
