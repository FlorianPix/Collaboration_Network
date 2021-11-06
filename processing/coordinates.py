from time import sleep
from random import randint
from typing import Optional
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderServiceError

from .model import Location, Coordinates


__user_agent = f"user_me_{randint(10000, 99999)}"
__geolocator = Nominatim(user_agent=__user_agent)
__SLEEP_SECONDS = 1


def get_coords(location: Location) -> Optional[Coordinates]:
    """get coords from location"""
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
        return get_coords(location)
    except GeocoderServiceError as error:
        print('CONNECTION REFUSED: GeocoderServiceError encountered.', error)
        return None


def get_location_coordinates(locations: list[Location]) -> dict[Location, Optional[Coordinates]]:
    """calculate dictionary containing coordinates of locations"""
    result: dict[Location, Coordinates] = {}
    for location in locations:
        if location not in result:
            result[location] = get_coords(location)
    return result
