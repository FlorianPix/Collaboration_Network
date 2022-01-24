"""fetching and dealing with coordinates"""
from random import randint
from time import sleep
from typing import Optional

from file_io import get_coords_pickle, save_coords_pickle
from geopy.exc import GeocoderServiceError, GeocoderTimedOut
from geopy.geocoders import Nominatim

from .model import Coordinates, Location, Paper
from .util import progressbar

__user_agent = f"user_me_{randint(10000, 99999)}"
__geolocator = Nominatim(user_agent=__user_agent)
__SLEEP_SECONDS = 1

# source: https://stackoverflow.com/a/60088688
def get_city_coords(location: Location) -> Optional[Coordinates]:
    """get coords of location utilizing the Nominatim api"""
    search = {
        'city': location.city,
        'country': location.country,
    }
    if location.state is not None:
        search['state'] = location.state
    try:
        result = __geolocator.geocode(search)
        return None if result is None else Coordinates(result.latitude, result.longitude)
    except GeocoderTimedOut:
        print('TIMED OUT: GeocoderTimedOut: Retrying...')
        sleep(randint(1*100, __SLEEP_SECONDS*100)/100)
        return get_city_coords(location)
    except GeocoderServiceError as error:
        print('CONNECTION REFUSED: GeocoderServiceError encountered.', error)
        return None

def calculate_paper_coordinates(papers: list[Paper]) -> dict[Location, Optional[Coordinates]]:
    """calculate dictionary containing coordinates of locations"""
    result: dict[Location, Optional[Coordinates]] = {}
    location_count = 0
    fail_count = 0
    try:
        known_coordinates = get_coords_pickle("coordinates.pkl")
    except FileNotFoundError:
        known_coordinates = {}
    for paper in progressbar(papers, "fetching coordinates: "):
        for location in paper.locations:
            location_count += 1
            if location not in result:
                if location in known_coordinates:
                    result[location] = known_coordinates[location]
                else:
                    tmp = get_city_coords(location)
                    if tmp is None:
                        fail_count += 1
                    result[location] = tmp
                    known_coordinates[location] = tmp
    if location_count > 0:
        print(f"Coords fail rate: {round(fail_count / location_count * 100, 2)}% (failed: {fail_count}, total: {location_count})")
    save_coords_pickle(known_coordinates, "coordinates.pkl")
    return result
