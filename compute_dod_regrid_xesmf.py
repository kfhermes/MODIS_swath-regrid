import os
import sys
from copy import deepcopy
from pathlib import Path

# import tqdm
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe


def set_value_to_nan(array, value):
    array = deepcopy(array)
    array[array==value] = np.nan
    return array

def get_files(path, suffix='.hdf', recursive=False):
    path = Path(path)
    if recursive:
        return list(path.rglob(f'*{suffix}'))
    else:
        return list(path.glob(f'*{suffix}'))


idir = sys.argv[1]
# idir = '/neodc/modis/data/MYD04_L2/collection61/2024/07/30'
odir = sys.argv[2]
# odir = 'directory/where/you/want/to/save/the/output'
os.makedirs(odir, exist_ok=True)
box = [-20, 50, 0, 40] # [lonmin, lonmax, latmin, latmax]


ifiles = get_files(path=idir, suffix='.hdf')
print(f"Found {len(ifiles)} files in {idir}")

# Search for files that intersect with the lonlatbox
num_points_in_box = np.zeros(len(ifiles))
# for idx in tqdm(range(len(ifiles))):
for idx in range(len(ifiles)):
    ds = xr.open_dataset(ifiles[idx], engine='netcdf4')

    swath_lon = ds['Longitude'].values
    swath_lat = ds['Latitude'].values

    lon_in_box = (swath_lon > box[0]) & (swath_lon < box[1])
    lat_in_box = (swath_lat > box[2]) & (swath_lat < box[3])

    lonlat_in_box = lon_in_box & lat_in_box
    num_points_in_box[idx] = np.sum(lonlat_in_box)

print("Found swath data within lonlatbox for indices:")
print(np.where(num_points_in_box>0)[0])


idx_data_avail = np.where(num_points_in_box>0)[0]
# for idx in tqdm(idx_data_avail):
for idx in idx_data_avail:
    print(f"Processing file {ifiles[idx]}")
    ds = xr.open_dataset(ifiles[idx], engine='netcdf4')
    aod = ds['Deep_Blue_Aerosol_Optical_Depth_550_Land'].values
    ae = ds['Deep_Blue_Angstrom_Exponent_Land'].values
    ssa = ds['Deep_Blue_Spectral_Single_Scattering_Albedo_Land'].values[1,:,:] # SSA at 470nm

    aod[ssa>=1.0] = np.nan # Set AOD to NaN where SSA is invalid

    dod = aod * (0.98 - 0.5089*ae + 0.0512*ae**2) # compute DOD

    # Get swath lon and lat
    swath_lon = ds['Longitude'].values
    swath_lat = ds['Latitude'].values

    # Define custom latitude and longitude grid
    target_lat = np.linspace(-90, 90, 721)
    target_lon = np.linspace(-180, 180, 1441)

    # Remove zero values by setting to NaN
    dod_zr = set_value_to_nan(dod, 0)

    # Create Array with 2D lat/lon coordinates
    grid_swath = {"lon": swath_lon, "lat": swath_lat}
    grid_target  = {"lat": target_lat, "lon": target_lon}

    # Create regridder
    regridder = xe.Regridder(ds_in=grid_swath, ds_out=grid_target, method="nearest_d2s")

    # Perform regridding
    dod_zr_rg = regridder(dod_zr)

    # To be used only with nearest_d2s: Compute number of mapped aod points
    dod_zr_rg_N = regridder((dod_zr>=0)*1)

    # Normalize accumulated points with number of mapped points
    dod_zr_rg = dod_zr_rg / dod_zr_rg_N

    # Remove zero values by setting to NaN
    dod_zr_rg = set_value_to_nan(dod_zr_rg, 0)
    
     # Save as NetCDF
    date_info = str(idir).split('/')
    YYYY = date_info[-3]
    mm = date_info[-2]
    dd = date_info[-1]

    data_info = str(ifiles[idx]).split('/')[-1].split('.')
    product = data_info[0]
    series = data_info[1]
    HHMM = data_info[2]
    identifier = data_info[4]

    patch_time = pd.to_datetime(f"{YYYY}-{mm}-{dd}T{HHMM[:2]}:{HHMM[2:]}:00")
    patch_time = np.array([patch_time])

    # Write netcdf output
    ds_out = xr.Dataset(
    {
        "DOD": (["time", "lat", "lon"], dod_zr_rg[None,:,:]),
    },
    coords={
        "time": patch_time,
        "lat": target_lat,
        "lon": target_lon,
    },
    )

    ds_out.attrs = {
        'title': 'Dust Optical Depth (DOD) derived from MODIS data',
        'summary': 'DOD calculated using MODIS AOD, SSA (470nm), and Angstrom exponent as described in Pu and Ginoux, 2018.',
        'references': 'https://doi.org/10.5194/acp-18-12491-2018',
    }

    encoding = {
    'DOD': {
        'dtype': 'float32', 
        'compression': 'zlib',
        'complevel': 4 ,
    }
    }
    
    fnm = f"{product}_{series}_{YYYY}{mm}{dd}_{HHMM}_{identifier}_DOD_0.25dgrid.nc"
    print(f"Saving file: {fnm}")

    os.makedirs(os.path.join(odir,f'{YYYY}{mm}'), exist_ok=True)
    
    ds_out.to_netcdf(os.path.join(odir,f'{YYYY}{mm}', fnm), encoding=encoding)