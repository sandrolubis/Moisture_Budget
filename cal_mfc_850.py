import xarray as xr
import numpy as np

def compute_mfc_centered(u, v, q):
    """
    Compute moisture flux convergence (MFC) at 850 hPa using centered finite difference.
    Assumes input arrays have shape (time, lat, lon) and regular 1D lat/lon.
    """
    Re = 6.371e6  # Earth radius in meters

    lat = u['lat']
    lon = u['lon']
    lat_rad = np.deg2rad(lat)
    dlat = np.deg2rad(lat[1] - lat[0])
    dlon = np.deg2rad(lon[1] - lon[0])

    # Compute moisture fluxes
    qu = q * u
    qv = q * v

    # Compute dx and dy (broadcast to shape)
    cos_lat = np.cos(lat_rad)
    dx = dlon * Re * cos_lat.broadcast_like(qu)
    dy = dlat * Re

    # Centered finite differences
    dqu_dx = (qu.roll(lon=-1) - qu.roll(lon=1)) / (2 * dx)
    dqv_dy = (qv.roll(lat=-1) - qv.roll(lat=1)) / (2 * dy)

    # Moisture flux convergence
    mfc = -(dqu_dx + dqv_dy)
    mfc.name = "mfc"
    mfc.attrs["long_name"] = "Moisture Flux Convergence"
    mfc.attrs["units"] = "1/s"
    return mfc

def main():
    # Load input data (plev=0 is the only level)
    u = xr.open_dataset("u850.nc")['var131'].isel(plev=0)
    v = xr.open_dataset("v850.nc")['var132'].isel(plev=0)
    q = xr.open_dataset("q850.nc")['var133'].isel(plev=0)

    # Compute MFC
    mfc = compute_mfc_centered(u, v, q)

    # Save to NetCDF
    mfc.to_netcdf("mfc_850.nc")
    print(" MFC successfully saved to mfc_850.nc")

if __name__ == "__main__":
    main()
