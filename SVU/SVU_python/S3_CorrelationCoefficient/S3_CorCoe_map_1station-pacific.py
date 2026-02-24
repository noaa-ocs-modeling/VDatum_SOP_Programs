#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import xarray as xr
import csv, glob, os
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# In[2]:


# ------------------------------
# Inputs (EDIT THESE PATHS/VALUES)
# ------------------------------
#MERGED_NC    = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/WaterDistance/S1_wdist_Ak/coef_full_mlw_matlabAligned_multi.nc"  # folder with shard files
MERGED_NC    = "/mnt/c/Users/mojgan.rostaminia/Documents/bin/Pacific_S3/coef_full_mtl_matlabAligned_multi.nc"
#GAGE_CSV = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/WaterDistance/S1_wdist_Ak/gageinfo1.csv"
GAGE_CSV ="/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/stations/Pacific_station_final/gageinfo1.csv"
#FORT14   = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/WaterDistance/S1_wdist_Ak/fort.14"
FORT14   ="/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/BaseMesh/ADCIRC_input/final_08062025/fort_geo.14"
DATUM    = "mlw"   # one of: mhhw, mhw, mlw, mllw, mtl, dtl

# choose how to select the station column in NC:
SELECT_BY = "station_id"           # "station_id" or "index"
PLOT_STATION_IDX   = 0             # used if SELECT_BY="index" (0 = first)
PLOT_STATION_ID_1B = None          # e.g., 9410230 (1-based ID as in CSV/NC); if None and SELECT_BY="station_id", uses first row in CSV

# zoom options
RADIUS = None        # e.g. 200_000.0, or None
XLIM   = None        # e.g. (x0, x1) or None
YLIM   = None        # e.g. (y0, y1) or None

# plotting knobs
MAX_POINTS = None    # e.g. 500000, or None for all
POINT_SIZE = 2.0


# In[3]:


# ================== HELPERS ==================
def pacific_wrap(lon, center=180.0):
    """Wrap longitudes so the Pacific is continuous (center near 180°)."""
    return ((lon - center + 180.0) % 360.0) - 180.0 + center
    
def read_fort14_nodes(path):
    with open(path, "r") as f:
        _ = f.readline()                # title
        parts = f.readline().split()
        if len(parts) < 2:
            raise RuntimeError("bad fort.14 header")
        a, b = int(parts[0]), int(parts[1])
        nnode = max(a, b)
        x = np.empty(nnode, dtype=float)
        y = np.empty(nnode, dtype=float)
        for i in range(nnode):
            sp = f.readline().split()
            if len(sp) < 4:
                raise RuntimeError(f"bad fort.14 node line {i+1}")
            x[i] = float(sp[1]); y[i] = float(sp[2])
    return x, y

def read_gage_csv(path):
    """Return arrays (station_ids_1b, node_ids_1b) from CSV (first two columns)."""
    sids, nids = [], []
    with open(path, "r", newline="") as f:
        rdr = csv.reader(f)
        _ = next(rdr, None)  # header
        for r in rdr:
            try:
                sids.append(int(r[0])); nids.append(int(r[1]))
            except Exception:
                pass
    return np.array(sids, dtype=np.int64), np.array(nids, dtype=np.int64)


# In[4]:


# ================== LOAD MERGED FILE ==================
ds = xr.open_dataset(MERGED_NC)
coef_da = ds["coef"]
print("coef dims:", coef_da.dims, "sizes:", {k:int(v) for k,v in coef_da.sizes.items()})

# Figure out dimension names (expects "station" and "node", any order). Fallback if unnamed.
dims = list(coef_da.dims)
has_named = ("station" in dims) and ("node" in dims)


# In[5]:


# Read gage CSV (we'll use it to place the star and to resolve station_id -> node_id)
csv_station_ids_1b, csv_node_ids_1b = read_gage_csv(GAGE_CSV)
if csv_station_ids_1b.size == 0:
    raise RuntimeError("gage CSV is empty or malformed")

# Decide which station to plot
if SELECT_BY.lower() == "station_id":
    sid_1b = PLOT_STATION_ID_1B if PLOT_STATION_ID_1B is not None else int(csv_station_ids_1b[0])
    # find station index in NC if station coordinate exists; otherwise infer from CSV order
    if has_named:
        st_coord = ds["station"].values
        hit = np.nonzero(st_coord == sid_1b)[0]
        if hit.size == 0:
            raise RuntimeError(f"station_id={sid_1b} not found in NC 'station' coord")
        station_idx = int(hit[0])
    else:
        # assume NC stations are stored in the same order the compute used (CSV order)
        # use position of station_id in CSV as index
        hit = np.nonzero(csv_station_ids_1b == sid_1b)[0]
        if hit.size == 0:
            raise RuntimeError(f"station_id={sid_1b} not found in CSV")
        station_idx = int(hit[0])
    print(f"[select] by station_id={sid_1b} → station_idx={station_idx}")
else:
    station_idx = int(PLOT_STATION_IDX)
    if has_named:
        sid_1b = int(ds["station"].values[station_idx])
    else:
        sid_1b = int(csv_station_ids_1b[station_idx])
    print(f"[select] by index={station_idx} → station_id={sid_1b}")


# In[6]:


with Dataset(MERGED_NC, "r") as nc:
    coef_vals = nc["coef"][:,0][:]            # shape (nnodes_total,)
#coef_vals = nc["coef"][:,43][:]   


# In[7]:


# Extract the 1-D coef array for that station over nodes
if has_named:
    cofst = coef_da.isel(station=station_idx)
    if cofst.ndim != 1:
        cofst = cofst.transpose("node")
    node_coord = ds.coords.get("node", xr.DataArray(np.arange(cofst.shape[0]), dims=["node"])).values
    coef_vals = cofst.values.astype(np.float32, copy=False)
else:
    # fallback if dims unnamed: assume one dim is stations ~ 301
    n0, n1 = coef_da.shape
    if min(n0, n1) < 2000:  # smaller dim is likely station (e.g., 301)
        if n0 < n1:
            coef_vals = coef_da.isel({coef_da.dims[0]: station_idx}).values.astype(np.float32, copy=False)
        else:
            coef_vals = coef_da.isel({coef_da.dims[1]: station_idx}).values.astype(np.float32, copy=False)
        node_coord = np.arange(coef_vals.shape[0], dtype=np.int64)
    else:
        raise RuntimeError("cannot infer station/node dims from merged file")
ds.close()

finite = np.isfinite(coef_vals)
print(f"[coef] length={coef_vals.size}, finite={finite.sum()}, min={np.nanmin(coef_vals):.3f}, max={np.nanmax(coef_vals):.3f}")


# In[8]:


# ================== MAP NODES TO COORDS ==================
x_all, y_all = read_fort14_nodes(FORT14)
# node_coord is expected to be 0-based absolute node indices.
xx = x_all[node_coord]
yy = y_all[node_coord]

# station star location: find matching row in CSV and get node_id_1b
hit = np.nonzero(csv_station_ids_1b == sid_1b)[0]
if hit.size == 0:
    # fallback: first row
    print(f"[warn] station_id {sid_1b} not found in CSV; using first row for star placement")
    station_node_1b = int(csv_node_ids_1b[0])
else:
    station_node_1b = int(csv_node_ids_1b[int(hit[0])])
# Save original X/Y for station, then apply wrap
sx_orig, sy_orig = x_all[station_node_1b - 1], y_all[station_node_1b - 1]

# --- APPLY THE PACIFIC WRAPPER HERE ---
xx = pacific_wrap(xx)
sx = pacific_wrap(sx_orig)
sy = sy_orig # Latitude (y) does not need wrapping
# --- END WRAPPER APPLICATION ---

print(f"[star] station_id={sid_1b}, node_id={station_node_1b}, at (x,y)=({sx:.1f},{sy:.1f})")


# In[9]:


# ================== ZOOM / FILTER ==================
mask = np.ones_like(coef_vals, dtype=bool)
if RADIUS is not None:
    # Note: sx and xx are now wrapped!
    dx, dy = xx - sx, yy - sy
    mask &= (dx*dx + dy*dy) <= (RADIUS*RADIUS)
if XLIM is not None:
    # Note: xx is now wrapped! You may need to adjust XLIM values (e.g., from [-180, -120] to [180, 240]) 
    # depending on how you define your desired window.
    mask &= (xx >= XLIM[0]) & (xx <= XLIM[1])
if YLIM is not None:
    mask &= (yy >= YLIM[0]) & (yy <= YLIM[1])

xx_plot, yy_plot, cc_plot = xx[mask], yy[mask], coef_vals[mask]

# ... (rest of the block remains the same)


# In[10]:


# ===== robust interactive + static plotting (20 distinct colors) =====
import numpy as np
import os, pathlib
from IPython.display import HTML, display, IFrame

# sanity: verify variables coming from your previous cells exist
need = ["xx", "yy", "coef_vals", "node_coord", "sx", "sy", "sid_1b"]
missing = [n for n in need if n not in globals()]
if missing:
    raise RuntimeError(f"Missing variables from earlier cells: {missing}")

# filter finite values
m = np.isfinite(xx) & np.isfinite(yy) & np.isfinite(coef_vals)
xxp, yyp, ccp, nodes_plot = xx[m], yy[m], coef_vals[m], node_coord[m]
print(f"[data] finite points: {xxp.size:,}")

# aggressive downsample for rendering (tune as needed)
MAX_POINTS = 300_000
if MAX_POINTS is not None and xxp.size > MAX_POINTS:
    rng = np.random.default_rng(0)
    sel = rng.choice(xxp.size, size=MAX_POINTS, replace=False)
    xxp, yyp, ccp, nodes_plot = xxp[sel], yyp[sel], ccp[sel], nodes_plot[sel]
    print(f"[data] downsampled to {MAX_POINTS:,} points")

# ---------- 1) write Plotly HTML and embed ----------
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib

# define 20 discrete bins from -1 to 1
bins = np.linspace(-1, 1, 21)  # 20 intervals
labels = 0.5 * (bins[:-1] + bins[1:])  # midpoints

# discretize coef values
inds = np.digitize(ccp, bins) - 1
inds = np.clip(inds, 0, 19)

# use matplotlib tab20 palette for 20 distinct colors
tab20 = matplotlib.cm.get_cmap("tab20", 20).colors
colors_hex = [matplotlib.colors.to_hex(c) for c in tab20]
colors_used = [colors_hex[i] for i in inds]

hover = [f"node={int(n)}<br>r={float(r):.2f}" for n, r in zip(nodes_plot, ccp)]

fig = go.Figure()
fig.add_trace(go.Scattergl(
    x=xxp, y=yyp,
    mode="markers",
    marker=dict(
        size=6,
        color=colors_used,
        showscale=False
    ),
    text=hover, hoverinfo="text",
    name="coef"
))
# station star
fig.add_trace(go.Scattergl(
    x=[sx], y=[sy],
    mode="markers",
    marker=dict(size=16, color="black", symbol="star"),
    name=f"Station {sid_1b}"
))
fig.update_layout(
    title=f"Interactive correlation map — station_id={sid_1b}",
    xaxis_title="X", yaxis_title="Y",
    dragmode="pan",
    template="plotly_white"
)
# add manual legend showing bin ranges
for i, mid in enumerate(labels):
    fig.add_trace(go.Scattergl(
        x=[None], y=[None],
        mode="markers",
        marker=dict(size=10, color=colors_hex[i]),
        legendgroup=f"bin{i}",
        showlegend=True,
        name=f"{bins[i]:.1f} to {bins[i+1]:.1f}"
    ))

html_path = "coef_interactive_map.html"
fig.write_html(html_path, include_plotlyjs="cdn")
print(f"[saved] {html_path} (open in a browser if it doesn't embed)")
try:
    display(IFrame(src=html_path, width="100%", height="700px"))
except Exception as e:
    print("[warn] iframe embed failed:", e)

# ---------- 2) also save a static PNG ----------
plt.figure(figsize=(14, 10))
sc = plt.scatter(xxp, yyp, c=inds, s=6, cmap="tab20", vmin=0, vmax=19)
if sx is not None:
    plt.scatter([sx], [sy], c="k", s=120, marker="*")
cbar = plt.colorbar(sc, ticks=np.linspace(0.5, 19.5, 20))
cbar.set_label("Correlation bins")
cbar.ax.set_yticklabels([f"{bins[i]:.1f}–{bins[i+1]:.1f}" for i in range(20)])
plt.xlabel("X"); plt.ylabel("Y")
plt.title(f"Static correlation map — station_id={sid_1b}")
png_path = "coef_static.png"
plt.savefig(png_path, dpi=180, bbox_inches="tight")
plt.close()
print(f"[saved] {png_path}")


# In[ ]:




