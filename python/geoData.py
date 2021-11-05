from json.decoder import JSONDecodeError
import time
import pandas as pd
from pandas.errors import EmptyDataError
from geopy.geocoders import Nominatim # OSM
from geopy.exc import GeocoderTimedOut
import json

def geolocate(cities):
    try:
        with open('data/knownCities.json', 'r') as file:
            knownCities = json.load(file)
    except JSONDecodeError as e:
        knownCities = {}

    coordinates = []
    geolocator = Nominatim(user_agent="johannalbrecht@protonmail.com", timeout=2)
    print("Geolocating cities...")    
    time0 = time.time()

    total = len(cities)
    i = 0
    unique = 0
    for alias in cities:
        i += 1
        if i % 1000 == 0:
            print("%d/%d done" %(i, total))
        if alias in knownCities.keys():
            coordinates.append(knownCities.get(alias))
        else:
            try:
                location = geolocator.geocode(alias)
                if location:
                    unique += 1
                    coordinates.append("(" + str(location.latitude) + ", " + str(location.longitude) + ")")
                    knownCities[alias] = "(" + str(location.latitude) + ", " + str(location.longitude) + ")"
            except GeocoderTimedOut as e:
                print("Error: geocode failed on input %s with message %s" %(alias, e))

    with open('data/cities.json', 'w') as file:
        json.dump(coordinates, file, indent=4)
    with open('data/knownCities.json', 'w') as file:
        json.dump(knownCities, file, indent=4)
    time1 = time.time()
    print("Done, %s seconds elapsed, %s unique cities out of %s total" %(time1 - time0, unique, total))    