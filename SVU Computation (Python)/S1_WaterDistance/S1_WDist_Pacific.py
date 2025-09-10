#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Developed by Lei Shi in matlab
Reformated by Liujuan Tang 01/29/2020 in matlab
Modified/Optimized and converted to python by Mojgan Rostaminia 09/10/2025
Ellipsoidal (Karney via pyproj) + Haversine options & Dijkstra python library

Create shortest water-path distance for each gauge.

USAGE EXAMPLES
--------------
  python S1_WDist_Pacific.py --gaugefile gageinfo.csv --meshfile fort.14
  python S1_WDist_Pacific.py -g gageinfo.csv -m fort.14 --iid 0 --outdir out_nc
  python S1_WDist_Pacific.py -g gageinfo.csv -m fort.14 --dry-thresh 0.0
  python S1_WDist_Pacific.py -g gageinfo.csv --list
  # NEW: choose distance engine (default: karney)
  python S1_WDist_Pacific.py -g gageinfo.csv -m fort.14 --distance haversine

INPUT
-----
gageinfo.csv  : text file with columns: station_id, node_id
fort.14

DEPENDENCIES
------------
Standard library: os, time, argparse, datetime
External: numpy, pandas, scipy, netCDF4, pyproj

Install with:
  conda install numpy pandas scipy netcdf4 pyproj 
"""

import os
import time
from datetime import datetime
import argparse
import numpy as np
import pandas as pd
import netCDF4 as nc
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

# ---- Ellipsoidal geodesics (Karney) via PROJ/pyproj ----
from pyproj import Geod
_GEOD = Geod(a=6378137.0, rf=298.257223563)  # explicit WGS-84


def log(msg: str):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def read_mesh_fort14(meshfile):
    """Read ADCIRC fort.14 (triangles only). Returns (title, numnod, numele, node_ids, x, y, z, tris_1based)."""
    t0 = time.time()
    with open(meshfile, "rt") as fid:
        title = fid.readline().strip()
        nums = np.loadtxt(fid, dtype=int, max_rows=1)
        numele, numnod = int(nums[0]), int(nums[1])

        BB = np.loadtxt(fid, dtype=float, max_rows=numnod)  # id, x, y, z
        node_ids = BB[:, 0].astype(np.int64)
        x = BB[:, 1].astype(np.float64)
        y = BB[:, 2].astype(np.float64)
        z = BB[:, 3].astype(np.float64)

        B = np.loadtxt(fid, dtype=int, max_rows=numele)     # id, 3, n1, n2, n3
        tris = B[:, [2, 3, 4]].astype(np.int64)             # keep 1-based indexing

    log(f"Read fort.14: {numnod:,} nodes, {numele:,} elements in {time.time()-t0:.1f}s")
    return title, numnod, numele, node_ids, x, y, z, tris


# ---------------- distance engines ----------------
def _wrap180(lon):
    """Wrap longitude to [-180, 180). Works with arrays."""
    return ((lon + 180.0) % 360.0) - 180.0


def _dist_km_haversine(lon1, lat1, lon2, lat2):
    """Spherical haversine (km) with R=6371 km. Seam-safe Δλ wrapping."""
    R = 6371.0
    rad = np.pi / 180.0
    dlon = _wrap180(lon2 - lon1) * rad
    dlat = (lat2 - lat1) * rad
    lat1r = lat1 * rad
    lat2r = lat2 * rad
    a = np.sin(dlat/2.0)**2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon/2.0)**2
    return 2.0 * R * np.arcsin(np.sqrt(a))


def _dist_km_karney_pyproj(lon1, lat1, lon2, lat2):
    """Ellipsoidal (WGS-84) via pyproj.Geod (Karney). Vectorized & fast. Returns kilometers."""
    _, _, s12 = _GEOD.inv(lon1, lat1, lon2, lat2)
    return s12 / 1000.0


def build_wet_graph(x, y, z, tris_1based, dry_thresh, distance_engine="karney"):
    """
    Build an undirected CSR graph of wet-only edges with geodesic weights (km).
    distance_engine: 'karney' (default) or 'haversine'
    Returns (G_csr, num_edges).
    """
    t0 = time.time()

    wet = z > dry_thresh
    T = tris_1based[
        wet[tris_1based[:, 0] - 1] &
        wet[tris_1based[:, 1] - 1] &
        wet[tris_1based[:, 2] - 1]
    ]
    if T.size == 0:
        raise ValueError("No wet triangles found with the given dry_thresh.")

    A = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]], axis=0)
    A = np.sort(A, axis=1)
    edges = np.unique(A, axis=0)
    i1 = edges[:, 0] - 1
    i2 = edges[:, 1] - 1

    lon1 = x[i1].astype(np.float64); lat1 = y[i1].astype(np.float64)
    lon2 = x[i2].astype(np.float64); lat2 = y[i2].astype(np.float64)

    if distance_engine == "karney":
        w_km = _dist_km_karney_pyproj(lon1, lat1, lon2, lat2)
        engine_name = "Karney (pyproj/PROJ, WGS-84)"
    elif distance_engine == "haversine":
        w_km = _dist_km_haversine(lon1, lat1, lon2, lat2)
        engine_name = "Haversine (sphere, R=6371 km)"
    else:
        raise ValueError(f"Unknown --distance option: {distance_engine}")

    row = np.concatenate([i1, i2])
    col = np.concatenate([i2, i1])
    dat = np.concatenate([w_km, w_km])
    G = coo_matrix((dat, (row, col)), shape=(x.size, x.size), dtype=np.float64).tocsr()

    log(f"Wet edges: {edges.shape[0]:,}. Graph built with {engine_name} in {time.time()-t0:.1f}s")
    return G, edges.shape[0]


def write_station_nc(outpath, station_id, node0, numnod, dist_km, predecessors):
    """Write one station's outputs to NetCDF."""
    dis = np.where(np.isfinite(dist_km), dist_km, 9.999e5).astype(np.float32)
    last_node = np.full(numnod, -1, dtype=np.int32)
    klayer = np.full(numnod, -1, dtype=np.int32)

    src = int(node0) - 1
    finite = np.isfinite(dist_km)
    last_node[finite] = predecessors[finite] + 1
    last_node[src] = 0

    order = np.argsort(dist_km)
    klayer[src] = 0
    for u in order:
        if u == src or not np.isfinite(dist_km[u]):
            continue
        p = predecessors[u]
        if p >= 0:
            klayer[u] = klayer[p] + 1

    with nc.Dataset(outpath, "w") as ds:
        ds.createDimension("mesh", numnod)
        ds.createDimension("node0", 1)
        v_dis = ds.createVariable("dis", "f4", ("mesh",))
        v_lst = ds.createVariable("last_node", "i4", ("mesh",))
        v_klr = ds.createVariable("klayer", "i4", ("mesh",))
        v_n0  = ds.createVariable("node0", "i4", ("node0",))

        v_dis[:] = dis
        v_lst[:] = last_node
        v_klr[:] = klayer
        v_n0[:]  = int(node0)

        v_dis.units = "km"
        v_dis.long_name = "shortest water-path distance from source gauge (over wet-mesh edges)"
        ds.title = "Shortest-path distance over wet ADCIRC mesh"
        ds.history = f"Created {datetime.now().isoformat()}"

    return outpath


def run_for_station(G, numnod, station_id, node0, outdir):
    src = int(node0) - 1
    if not (0 <= src < numnod):
        raise IndexError(f"node0 {node0} out of range for mesh with {numnod} nodes")

    t1 = time.time()
    dist_km, predecessors = dijkstra(G, directed=False, indices=src, return_predecessors=True)
    fn = os.path.join(outdir, f"station{int(station_id)}.nc")
    write_station_nc(fn, int(station_id), int(node0), numnod, dist_km, predecessors)
    return fn, time.time() - t1


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute shortest water-path distances over an ADCIRC mesh (fort.14) for one or more gauges."
    )
    p.add_argument("--gaugefile", "-g", default="gageinfo.csv",
                   help="CSV with columns: station_id,node_id (default: gageinfo.csv)")
    p.add_argument("--meshfile", "-m", default="fort.14",
                   help="ADCIRC mesh file (default: fort.14)")
    p.add_argument("--dry-thresh", type=float, default=0.1,
                   help="Nodes with z <= dry_thresh are dry (default: 0.1)")
    p.add_argument("--iid", type=int,
                   help="If provided, only process row index iid (0-based) from gaugefile")
    p.add_argument("--outdir", default=".",
                   help="Directory for outputs (default: current directory)")
    p.add_argument("--list", action="store_true",
                   help="List gauges parsed from gaugefile and exit")
    p.add_argument("--distance", choices=["karney", "haversine"], default="karney",
                   help="Edge-length metric: 'karney' (ellipsoidal WGS-84, default) or 'haversine' (sphere R=6371).")
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()
    log("Starting")

    df = pd.read_csv(args.gaugefile)
    if not {"station_id", "node_id"}.issubset(df.columns):
        raise ValueError("gaugefile must have columns: station_id,node_id")

    if args.list:
        print(df.head(20).to_string(index=True))
        print(f"... ({len(df)} total rows)")
        return

    if args.iid is not None:
        if args.iid < 0 or args.iid >= len(df):
            raise IndexError(f"iid={args.iid} out of range; file has {len(df)} rows")
        df = df.iloc[[args.iid]].reset_index(drop=True)

    station_ids = df["station_id"].to_numpy(np.int64)
    node0s      = df["node_id"].to_numpy(np.int64)
    log(f"Loaded {len(df)} gauge(s) from {args.gaugefile}")

    title, numnod, numele, node_ids, x, y, z, tris = read_mesh_fort14(args.meshfile)
    log(f"Mesh title: {title}")
    log(f"x range: {x.min():.6f} .. {x.max():.6f} | y range: {y.min():.6f} .. {y.max():.6f} | z range: {z.min():.3f} .. {z.max():.3f}")

    G, n_edges = build_wet_graph(x, y, z, tris, args.dry_thresh, distance_engine=args.distance)

    os.makedirs(args.outdir, exist_ok=True)
    ok = 0
    for k, (sid, node0) in enumerate(zip(station_ids, node0s), 1):
        try:
            fn, elapsed = run_for_station(G, numnod, sid, node0, args.outdir)
            ok += 1
            log(f"[{k:02d}/{len(df)}] station {sid} → {fn} | {elapsed:.1f}s")
        except Exception as e:
            log(f"[{k:02d}/{len(df)}] station {sid} FAILED: {e}")

    log(f"Done. Wrote {ok}/{len(df)} station files. Total time {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
