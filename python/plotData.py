import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
from pyproj import CRS
import json

def plotCities():
    print("Plotting...")
    with open("data/cities.json", "r") as file:
        coordinates = json.load(file)
    df = pd.DataFrame.from_dict(coordinates)
    geometry = [Point(float(x[1:-1].split(", ")[1]), float(x[1:-1].split(", ")[0])) for x in df[0]]
    crs = CRS('EPSG:4326') # coordinate system
    geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    countries_map = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    f, ax = plt.subplots(figsize=(16, 16))
    countries_map.plot(ax=ax, alpha=0.4, color='grey')
    geo_df['geometry'].plot(ax=ax, markersize=30, color='b', marker='^', alpha=.2)
    plt.show()