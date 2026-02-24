#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import warnings

warnings.filterwarnings("ignore")

# --------------------------------------------------------------
# PHASE 2 NODE-CHUNK WORKER
# --------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runid", required=True)
    ap.add_argument("--datum", required=True,
                    choices=["mllw", "mhhw", "mhw", "mlw", "mtl", "dtl"])
    ap.add_argument("--chunk-id", required=True, type=int)
    ap.add_argument("--nchunks", required=True, type=int)
    ap.add_argument("--path_pre", required=True)
    args = ap.parse_args()

    runid = args.runid
    datum = args.datum
    kid = args.chunk_id
    nchunks = args.nchunks
    path_pre = args.path_pre

    # Fixed paths
    indata = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/indata/"
    wdist_dir = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S1_Wdist/"
    kcoef_dir = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S3_CC/Pacific/"

    # -----------------------------
    # 1. Load dm_<datum>.nc
    # -----------------------------
    dm_file = os.path.join(path_pre, f"dm_{datum}.nc")
    try:
        with Dataset(dm_file, 'r') as nc:
            dm = nc.variables['z'][:].flatten()
    except Exception as e:
        print(f"Error loading {dm_file}: {e}")
        return

    n = len(dm)

    # -----------------------------
    # 2. Load Observation CSV
    # -----------------------------
    obs_file = "/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/Obs_Pacific_SVU.csv"
    print(f"Loading observations from {obs_file}...")
    
    try:
        df = pd.read_csv(obs_file)
        df.columns = df.columns.str.lower()
        
        # --- REMOVED SLICING: Reads entire CSV now ---

        required = ['id', 'node', 'rms', datum]
        if not all(col in df.columns for col in required):
            print(f"Error: CSV missing columns. Found: {df.columns.tolist()}")
            return

        sid_nodes = df['node'].values.astype(int)
        obs_ids = df['id'].values.astype(int)
        obs_vals = df[datum].values
        obs_rms = df['rms'].values

    except Exception as e:
        print(f"Error loading CSV {obs_file}: {e}")
        return

    # -----------------------------
    # Build PMOE List
    # -----------------------------
    pmoe = []
    for i in range(len(df)):
        val = obs_vals[i]
        r = obs_rms[i]
        node_idx = sid_nodes[i]
        station_id = obs_ids[i]

        if np.isnan(val) or val == 0: continue
        if np.isnan(r) or r > 999 or r < -9: continue
            
        pmoe.append((station_id, node_idx, float(val), float(r)))

    m = len(pmoe)
    print(f"Processed {m} valid observations.")

    # Convert to arrays
    bcpsi = np.array([p[1] for p in pmoe])
    ostn = np.array([p[2] for p in pmoe]).reshape(m, 1)
    uncer = np.array([p[3] for p in pmoe]).reshape(m, 1)
    
    if np.max(bcpsi - 1) >= n:
        print(f"ERROR: Max station index ({np.max(bcpsi)}) exceeds grid size ({n}).")
        return

    mstn = dm[bcpsi - 1].reshape(m, 1)
    estn = ostn - mstn

    o_error = float(np.sqrt(np.mean(uncer**2)))
    m_error = float(np.std(estn, ddof=1))
    Lr = 2.0

    # ----------------------------------------
    # Chunking
    # ----------------------------------------
    ns = n // nchunks
    i0 = kid * ns
    i1 = n if kid == nchunks - 1 else (kid + 1) * ns

    dm_s = dm[i0:i1]
    nc_s = i1 - i0

    # ----------------------------------------
    # LOAD COEFFICIENTS
    # ----------------------------------------
    kcoef_file = os.path.join(kcoef_dir, f"coef_full_{datum}_matlabAligned_multi.nc")
    
    kc_chunk_all = np.ones((nc_s, 1), dtype='f8') 
    coef_station_map = {} 
    
    print(f"Reading coefficient chunk from {kcoef_file}...")
    try:
        with Dataset(kcoef_file, 'r') as nc_kc:
            # 1. Read coef chunk (rows=current nodes, cols=all stations in file)
            kc_chunk_all = nc_kc.variables['coef'][i0:i1, :]
            
            # 2. Map Station IDs to Column Index
            if 'station' in nc_kc.variables:
                file_stations = nc_kc.variables['station'][:].astype(int)
                # Map: StationID -> ColumnIndex
                coef_station_map = {sid: idx for idx, sid in enumerate(file_stations)}
            else:
                print("Warning: 'station' variable missing in coef file.")

    except Exception as e:
        print(f"Warning: Could not load coef ({e}). Using ones.")
        # Fallback to avoid crash loop, using m (number of processed obs)
        kc_chunk_all = np.ones((nc_s, m), dtype='f8')

    # ----------------------------------------
    # Output file (NetCDF)
    # ----------------------------------------
    outfile = os.path.join(path_pre, f"SVU_partial_{datum}_{runid}_chunk{kid:04d}.nc")

    with Dataset(outfile, "w", format="NETCDF4") as nc:
        nc.createDimension("node", nc_s)
        nc.createDimension("station", m)
        nc.createDimension("one", 1)

        nc.createVariable("dm", "f8", ("node",))[:] = dm_s
        nc.createVariable("bcpsi", "i4", ("station",))[:] = bcpsi
        nc.createVariable("mstn", "f8", ("station",))[:] = mstn[:, 0]
        nc.createVariable("ostn", "f8", ("station",))[:] = ostn[:, 0]
        nc.createVariable("estn", "f8", ("station",))[:] = estn[:, 0]
        nc.createVariable("uncer", "f8", ("station",))[:] = uncer[:, 0]
        nc.createVariable("o_error", "f8")[:] = o_error
        nc.createVariable("m_error", "f8")[:] = m_error
        nc.createVariable("Lr", "f8")[:] = Lr
        
        var = nc.createVariable("varname", str, ("one",))
        var[0] = datum

        covar = nc.createVariable("covar0", "f8", ("node", "station"))

        # ----------------------------------------
        # Compute Covariance Columns
        # ----------------------------------------
        for j, (sid, snode, _, _) in enumerate(pmoe):

            # 1. Load Distance Chunk
            dist_file = os.path.join(wdist_dir, f"station{sid}.nc")
            md = np.full(nc_s, np.inf)
            
            try:
                with Dataset(dist_file, 'r') as nc_dist:
                    if 'dis' in nc_dist.variables:
                        md = nc_dist.variables['dis'][i0:i1] / 111.0
            except Exception:
                pass 

            # 2. Match Station ID (CSV) to Station ID (NetCDF)
            if sid in coef_station_map:
                col_idx = coef_station_map[sid]
                # Ensure the mapped column index is valid for the loaded chunk
                if col_idx < kc_chunk_all.shape[1]:
                    kc = kc_chunk_all[:, col_idx]
                else:
                    kc = np.ones(nc_s)
            else:
                # Station from CSV not found in NetCDF file -> No correlation adjustment
                kc = np.ones(nc_s)

            # 3. Calculations
            d_ref = dm[snode - 1]
            if d_ref == 0: d_ref = 1e-4

            ratio = np.abs(dm_s / d_ref)
            mask_r = ratio > 1
            ratio[mask_r] = 1.0 / ratio[mask_r]

            cov = (m_error**2) * np.exp(-md / Lr) * (kc * ratio)

            covar[:, j] = cov

    print(f"Chunk {kid} complete: {outfile}")

if __name__ == "__main__":
    main()
