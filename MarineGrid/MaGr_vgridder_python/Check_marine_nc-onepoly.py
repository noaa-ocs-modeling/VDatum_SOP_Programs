#!/usr/bin/env python
# coding: utf-8

import netCDF4 as nc
import numpy as np

#nc_path = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/marine_PA_W_CNMI_01.nc"
nc_path = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/marine_PA_C_Ocean_01.nc"
ds = nc.Dataset(nc_path, "r")
field = ds.variables['field'][:]
unique_vals = np.unique(field)

print(f"Unique values in {nc_path}: {unique_vals}")
ds.close()





