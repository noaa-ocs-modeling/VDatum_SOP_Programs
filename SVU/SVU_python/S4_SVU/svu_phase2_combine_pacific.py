#!/usr/bin/env python3
import os
import argparse
import numpy as np
from netCDF4 import Dataset

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runid", required=True)
    ap.add_argument("--datum", required=True)
    ap.add_argument("--nchunks", required=True, type=int)
    ap.add_argument("--path_pre", required=True)
    args = ap.parse_args()

    runid = args.runid
    datum = args.datum
    nchunks = args.nchunks
    path_pre = args.path_pre

    f0 = os.path.join(path_pre,
        f"SVU_partial_{datum}_{runid}_chunk0000.nc")

    nc0 = Dataset(f0)
    m = len(nc0.dimensions["station"])
    dm_all = []
    cov_all = []
    bcpsi = nc0["bcpsi"][:]
    mstn = nc0["mstn"][:]
    ostn = nc0["ostn"][:]
    estn = nc0["estn"][:]
    uncer = nc0["uncer"][:]
    o_error = float(nc0["o_error"][:])
    m_error = float(nc0["m_error"][:])
    Lr = float(nc0["Lr"][:])
    nc0.close()

    for kid in range(nchunks):
        f = os.path.join(path_pre,
            f"SVU_partial_{datum}_{runid}_chunk{kid:04d}.nc")
        nc = Dataset(f)
        dm_all.append(nc["dm"][:])
        cov_all.append(nc["covar0"][:])
        nc.close()

    dm = np.concatenate(dm_all)
    covar = np.concatenate(cov_all)

    outfile = os.path.join(path_pre,
        f"SVU_input_{datum}_{runid}.nc")

    if os.path.exists(outfile):
       os.remove(outfile)
    nc = Dataset(outfile, "w", format="NETCDF4")
    n = len(dm)

    nc.createDimension("node", n)
    nc.createDimension("station", m)

    nc.createVariable("dm", "f8", ("node",))[:] = dm
    nc.createVariable("bcpsi", "i4", ("station",))[:] = bcpsi
    nc.createVariable("mstn", "f8", ("station",))[:] = mstn
    nc.createVariable("ostn", "f8", ("station",))[:] = ostn
    nc.createVariable("estn", "f8", ("station",))[:] = estn
    nc.createVariable("uncer", "f8", ("station",))[:] = uncer

    nc.createVariable("o_error", "f8")[:] = o_error
    nc.createVariable("m_error", "f8")[:] = m_error
    nc.createVariable("Lr", "f8")[:] = Lr

    nc.createDimension("one", 1)
    var = nc.createVariable("varname", str, ("one",))
    var[0] = datum

    nc.createVariable("covar0", "f8", ("node", "station"))[:] = covar

    nc.close()
    print(f"COMBINED → {outfile}")


if __name__ == "__main__":
    main()
