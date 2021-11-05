import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point

def plotCities():
    print("Plotting...")
    df = pd.read_csv("cities.csv")
    geometry = [Point(float(x[1:-1].split(", ")[1]), float(x[1:-1].split(", ")[0])) for x in df['Coordinates']]
    crs = {'init': 'epsg:4326'} # coordinate system
    geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    countries_map = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    f, ax = plt.subplots(figsize=(16, 16))
    countries_map.plot(ax=ax, alpha=0.4, color='grey')
    geo_df['geometry'].plot(ax=ax, markersize=30, color='b', marker='^', alpha=.2)
    plt.show()