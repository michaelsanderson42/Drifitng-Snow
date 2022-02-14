"""
Estimates the number of days in the future which could have drifting snow
and blizzards, from the UKCP18 Local projections.
"""

import os
import getpass
import numpy as np
import iris
from iris.time import PartialDateTime
from iris.analysis.cartography import rotate_pole
import matplotlib.pyplot as plt
from whichbox import whichbox


def get_ws_in_ms(ws, units='kn'):
    """
    Converts wind speed in knots or mph to metres per second.
    1 knot (kn, kt or kts) = 1.15077945 miles per hour (mph) = 1.852 kilometer per hour (kph)
      = 1.68780986 foot / second (ft/s) = 0.514444444 meters per second (m/s).
    1 nautical mile is now defined as exactly 1.852 km.
    """

    if units in ['kn', 'kt', 'kts', 'knot', 'knots']:
    # Convert wind speed in knots to metres per second
        wspeed = ws * 0.514444444
    elif units == 'mph':
    # Convert wind speed in miles per hour to metres per second
        wspeed = ws * (1760 * 3 * 12 * 2.54) / (3600 * 100)

    return wspeed


def make_ukcp_local_datadir(ensemble_member, var_name):
    """
    Constructs the name of the directory containing the UKCP Local data

    ensemble_member -- The UKCP ensemble member, an integer between 1 and 15.
    var_name -- The CF-name of the UKCP variable
    """

    if var_name == 'wsgmax10m':
        datadir = '/scratch/hadmi/UKCP18/ukcp18_data/'
    else:
        datadir = f"/project/ukcp/land-cpm/uk/2.2km/rcp85/{ensemble_member:02d}/{var_name}/day/v20210615/"

    return datadir


def make_ukcp_filename(ensemble_member, var_name, the_year):
    """
    Constructs the UKCP Local filename for the given ensemble member and variable,
    including the directory. The daily data are stored in 1 year chunks.

    ensemble_member -- The UKCP ensemble member, an integer between 1 and 15.
    var_name -- The CF-name of the UKCP variable
    the_year -- The year required.
    """

    datadir = make_ukcp_local_datadir(ensemble_member, var_name)
    fhead = f"{var_name}_rcp85_land-cpm_uk_2.2km_{ensemble_member:02d}_day"
    fend = "{:4d}1201-{:4d}1130.nc".format(the_year, the_year+1)

    return os.path.join(datadir, '_'.join([fhead, fend]))


def ukcp_callback(cube, field, filename):

# Remove 2D coordinates and other coords not needed.
    cube.attributes = {}
    cube.remove_coord('latitude')
    cube.remove_coord('longitude')
    cube.remove_coord('ensemble_member_id')
    cube.remove_coord('yyyymmdd')


def read_ukcp18_local(ensemble_member, var_name, the_year, lon_limits, lat_limits):
    """
    Reads in UKCP18 Local data for the extended winter season (October to May)

    ensemble_member -- The UKCP ensemble member, an integer between 1 and 15.
    var_name -- The CF-name of the UKCP variable
    the_year -- The year required. Data for October [year-1] to May [year] are loaded.
    lon_limits -- Only load data for grid boxes within this limit
    lat_limits -- Only load data for grid boxes within this limit
    """

    filename_y1 = make_ukcp_filename(ensemble_member, var_name, the_year-1)
    filename_y2 = make_ukcp_filename(ensemble_member, var_name, the_year)
    print(the_year)
    print(lon_limits)
    print(lat_limits)

# Constraints for the given lon/lat limits
    lon_con = iris.Constraint(grid_longitude = lambda l: lon_limits[0] < l.point < lon_limits[1])
    lat_con = iris.Constraint(grid_latitude = lambda l: lat_limits[0] < l.point < lat_limits[1])
    e_con = iris.Constraint(ensemble_member = ensemble_member)
    cube_1 = iris.load_cube(filename_y1, constraint=e_con & lat_con & lon_con, callback=ukcp_callback)
    cube_2 = iris.load_cube(filename_y2, constraint=e_con & lat_con & lon_con, callback=ukcp_callback)
    cube = iris.cube.CubeList([cube_1, cube_2]).concatenate_cube()

# Now constrain on time
    pdt1 = PartialDateTime(year=the_year, month=10, day=1)
    pdt2 = PartialDateTime(year=the_year+1, month=5, day=30)

    return cube.extract(iris.Constraint(time = lambda t: pdt1 <= t.point <= pdt2))


def make_empty_cube(nlats, nlons, var_name, the_year, lat_coord, lon_coord):
    """
    Makes an empty cube (i.e. a cube whose data are all zeros)
    """

    tcoord = iris.coords.AuxCoord(the_year, long_name='year', var_name='year', units='1')
    data = np.zeros((nlats, nlons), dtype=np.int)

    cube = iris.cube.Cube(data, long_name = f'days with {var_name} above/below threshold',
        var_name=var_name, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)])
    cube.add_aux_coord(tcoord)

    return cube


def find_snow_drift_days(the_year, tasmin, gust, snw):
    """
    Identifies days on which drifting snow or blizzards could occur, using thresholds
    from the forecasters' handbook (3rd Edn).
    https://digital.nmla.metoffice.gov.uk/IO_74afe7b9-24f7-4de9-96dd-20f4d5ac4094/
    """

# Thresholds for drifting snow, from the 3rd Edn of the Forecaster's Handbook.
# SWE of 0.02 mm used in the revised UKCP Local report, by Kendon et al. (2021), but is too low.
# 10 mm from critical value of Roesch et al. 2001, Clim Dyn.
    tlo = 0.0        #  Temperature must be below this threshold for drifting snow
    gmin_kn = 12.0   #  Wind speed in knots above which drifting snow occurs.
    swe_min = 10.0   #  Snow water equivalent (in mm) must be greater than this limit for there to be sufficient lying snow

# Constants for estimating snow cover fraction from snow water equivalent. SWE must be in metres.
# See Roesch et al. 2001, Clim Dyn. and Yang et al. (1997), J Climate, pp.353-373
    a = 100.0
    b = 0.95

# Convert speed in knots to m/s.
    gmin = get_ws_in_ms(gmin_kn, units='kn')

    ntimes, nlats, nlons = tasmin.shape

# Get coordinate for output cubes
    lon_coord = tasmin.coord('grid_longitude')
    lat_coord = tasmin.coord('grid_latitude')

# Set up output cubes for number of days with tasmin < tlo, wind gusts greater than gmin, number of days with lying
# snow, number of days with a snow fraction > 0.5, and number of days with drifting snow.
    d_below = make_empty_cube(nlats, nlons, tasmin.var_name, the_year, lat_coord, lon_coord)
    d_windy = make_empty_cube(nlats, nlons, gust.var_name, the_year, lat_coord, lon_coord)
    d_lying_snow = make_empty_cube(nlats, nlons, snw.var_name, the_year, lat_coord, lon_coord)
    d_snow_frac = make_empty_cube(nlats, nlons, 'snow_frac', the_year, lat_coord, lon_coord)
    d_drift = make_empty_cube(nlats, nlons, 'drift', the_year, lat_coord, lon_coord)

    for j in list(range(nlats)):
        for i in list(range(nlons)):
            # Find all days on which each criteria for drifting snow is met, and days when all three are met simultaneously
            wt = (tasmin.data[:,j,i] < tlo)
            wg = (gust.data[:,j,i] > gmin)
            ws = (snw.data[:,j,i] > swe_min)
            wd = (wt & wg & ws)
            snow_fracs = b * np.tanh(a * snw.data[:,j,i]/1000)
            wf = (snow_fracs >= 0.5)  #  Obs of snow cover require at least 50% coverage.

            d_below.data[j,i] = sum(wt)
            d_windy.data[j,i] = sum(wg)
            d_lying_snow.data[j,i] = sum(ws)
            d_snow_frac.data[j,i] = sum(wf)
            d_drift.data[j,i] = sum(wd)

    return iris.cube.CubeList([d_below, d_windy, d_lying_snow, d_snow_frac, d_drift])
    

def count_snow_drift_days(ensemble_member, year_range, lon_limits, lat_limits):

    dout = '/scratch/hadmi/UKCP18/processed/'

    # UKCP18 variable names
    var_name_t = 'tasmin'
    var_name_g = 'wsgmax10m'
    var_name_s = 'snw'

# Analyse data for October to May, so start at year+1
    for the_year in list(range(year_range[0]+1, year_range[1])):
        print('Loading data for year', the_year)

        cube_t = read_ukcp18_local(ensemble_member, var_name_t, the_year, lon_limits, lat_limits)
        cube_s = read_ukcp18_local(ensemble_member, var_name_s, the_year, lon_limits, lat_limits)
        print(np.min(cube_s.data), np.max(cube_s.data))

        cube_g10m = read_ukcp18_local(ensemble_member, var_name_g, the_year, lon_limits, lat_limits)
        # Estimate 2m wind gusts from 10m gusts
        cube_g = calculate_2m_gusts(cube_g10m)

        clist = find_snow_drift_days(the_year, cube_t, cube_g, cube_s)

        fout = os.path.join(dout, f'drift_rcp85_land-cpm_uk_2.2km_{ensemble_member:02d}_{the_year:4d}.nc')
        iris.save(clist, fout)


def calculate_2m_gusts(cube):
    """
    Calculates gust speeds at 2m above the surface, assuming log wind profile is valid.

    cube -- Iris cube containing 10m wind gusts
    """

    z0 = 0.03  # Roughness length (m) for open moorland, from Hammond et al. (2011), https://doi.org/10.1002/met.273
    heather_height = 0.30   # Heather on moorland assumed to be about 30 cm high
    z2 = 2     # Height of weather station
    z1 = 10    # Height of modelled wind gust speed

    # Zero plane displacement, m.
    d = heather_height * 0.67

# Calculate wind speed at 2m from speed at 10 m using:
# u2 = u10 * ln((z2 - d) / z0) / ln((z1 - d) / z0)
# where z2 = 2 m, z1 = 10 m
    
    fac = np.log((z2 - d) / z0) / np.log((z1 - d) / z0)

    cube_out = cube.copy()
    cube_out.data *= fac
    cube_out.var_name = 'wsgmax2m'

    return cube_out


def get_area_limits(ensemble_member, the_lat, the_lon, the_year, radius):

# Read in sample data on the rotated grid
    var_name = 'tasmin'
    filename = make_ukcp_filename(ensemble_member, var_name, the_year)
    e_con = iris.Constraint(ensemble_member = ensemble_member)
    cube = iris.load_cube(filename, constraint=e_con, callback=ukcp_callback)

# Get the data size, then the grid spacing
    ntimes, nlats, nlons = cube.shape
    nx = nlons // 2
    ny = nlats // 2
    grid_latitudes = cube.coord('grid_latitude').points
    grid_longitudes = cube.coord('grid_longitude').points
    dx, = np.diff(grid_longitudes[nx:nx+2])
    dy, = np.diff(grid_latitudes[ny:ny+2])

# Convert the coordinates of the snow shed to rotated coordinates
    pole_lon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    pole_lat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    lons = np.array([the_lon])
    lats = np.array([the_lat])
    rotated_lons, rotated_lats = rotate_pole(lons, lats, pole_lon, pole_lat)
    rotated_lons += 360

# Find the indices of grid box containing the point of interest
    i_pts = whichbox(cube, [rotated_lons[0], rotated_lats[0]],
        xcoord_name='grid_longitude', ycoord_name='grid_latitude')

# Coordinates of centre of grid box containing the point of interest
    clon = grid_longitudes[i_pts[0]]
    clat = grid_latitudes[i_pts[1]]

# Get edges of an area of n x n boxes around this point, where n = 2*radius + 1
    lon_l = clon - dx * (radius + 0.1)
    lon_r = clon + dx * (radius + 0.1)
    lat_b = clat - dy * (radius + 0.1)
    lat_t = clat + dy * (radius + 0.1)

    return (lon_l, lon_r), (lat_b, lat_t)


if __name__ == '__main__':

# Coordinates of the Cruach Snow Shed (from the customer).
    clat = 56.70308
    clon = -4.58581

# Radius in grid boxes of an area around the snow shed to study
    radius = 2

# Get the limits of a box around the shed. radius = 2 means a 5 x 5 grid of boxes
    ensemble_member = 1
    year = 1981
    lon_limits, lat_limits = get_area_limits(ensemble_member, clat, clon, year, radius)

#   year_ranges = [[1980, 2000], [2020, 2040], [2060, 2080]]  #  Models run from December to November
    year_ranges = [[2020, 2040], [2060, 2080]]  #  Models run from December to November

    for year_range in year_ranges:
        count_snow_drift_days(ensemble_member, year_range, lon_limits, lat_limits)

