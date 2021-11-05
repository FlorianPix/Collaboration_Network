import time
import pandas as pd
from geopy.geocoders import Nominatim # OSM
from geopy.exc import GeocoderTimedOut

def geolocate(cities):
    geolocator = Nominatim(user_agent="johannalbrecht@protonmail.com", timeout=2)
    locations = []
    knownLocations = {}
    print("Locating cities...")
    time0 = time.time()
    for city in cities:
        if city in knownLocations.keys():
            locations.append(knownLocations.get(city))
        else:
            try:
                location = geolocator.geocode(city)
                if location:
                    print(location)
                    locations.append(location)
                    knownLocations[city] = location
            except GeocoderTimedOut as e:
                print("Error: geocode failed on input %s with message %s" %(city, e))
    pd.DataFrame(locations, columns=['City Name', 'Coordinates']).to_csv("cities.csv")
    time1 = time.time()
    print("Done, %s seconds elapsed" %(time1 - time0))