from time import sleep
from random import randint
from typing import Optional
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderServiceError

from .model import Location, Coordinates, Paper


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
    index = 0
    for paper in papers:
        print(index)
        index += 1
        for location in paper.locations:
            if location not in result:
                result[location] = get_city_coords(location)
    return result
