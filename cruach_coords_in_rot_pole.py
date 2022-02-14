import iris
from iris.analysis.cartography import rotate_pole
import numpy as np

clat = 56.70308
clon = -4.58581

lons = np.array([clon])
lats = np.array([clat])

pole_lon = 177.5
pole_lat = 37.5

rotated_lons, rotated_lats = rotate_pole(lons, lats, pole_lon, pole_lat)

print(rotated_lats, 360.0+rotated_lons)
