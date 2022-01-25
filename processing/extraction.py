"""Extract location data out of string"""
from typing import Any, Optional, Tuple
import json
import geograpy
import geotext

from .model import Location, Paper

__countries = json.load(open("data/dictionaries/countries.json", 'rb'))

def __extract_state(target: str) -> Optional[Tuple[str, str]]: # (state, country)
    for country, data in __countries.items():
        states = data.get("states", {})
        for state in states:
            if state.startswith(target) or target.startswith(state):
                return (state, country)
        for state, alternatives in states.items():
            for alternative in alternatives:
                if target.startswith(alternative) or target.endswith(alternative):
                    return (state, country)
    return None


def __extract_country(target: str) -> Optional[Tuple[Optional[str], str]]: # (state, country)
    for country in __countries:
        if country.startswith(target) or target.startswith(country):
            return None, country
    for country, data in __countries.items():
        for alt in data["alt"]:
            if alt.startswith(target) or target.startswith(alt):
                return None, country
    for country, data in __countries.items():
        code = data["code"]
        if target.startswith(code) or target.endswith(code):
            return None, country
    result = __extract_state(target)
    if result is not None:
        return result
    # print(f"Unknown country: '{target}'")
    return None

def __clean_city(target: str) -> Optional[str]:
    parts = target.split(" ")
    result = []
    for part in parts:
        if not any(char.isdigit() or char == '@' for char in part):
            result.append(part.strip())
    result = " ".join(result)
    if not result:
        return None
    return result

# there's a lot of room for optimization here
def get_location_naive(text: str) -> Optional[Location]:
    """extract location using naive parsing"""
    if text.endswith('.'):
        text = text[:-2]
    parts = text.split(',')
    if len(parts) < 3:
        parts = text.split(' ')
    if len(parts) < 3:
        # print(f"Bad format: {text}")
        return None
    target = parts[-1].strip().lower()
    result = __extract_country(target)
    if result is None:
        return None
    state, country = result
    city = parts[-2].strip().lower()
    # in case the second last position contains a state instead of city
    result = __extract_state(city)
    if result is not None:
        state2, country2 = result
        if country2 == country:
            city = parts[-3].strip().lower()
            state = state2
    city = __clean_city(city)
    if city is None:
        return None
    return Location(city, state, country)

# this is extremely slow and does not work very well either;
# need to run `geograpy-nltk` in order for it to work;
# if it still does not work,
# save locations.db database from https://github.com/somnathrakshit/geograpy3/wiki
# in ~/.geograpy3
def get_location_geograpy(text: str) -> Optional[Location]:
    """extract location using geograpy"""
    places = geograpy.get_geoPlace_context(text=text)
    if not places.cities or not places.countries:
        return None
    region = places.regions[0] if places.regions else None
    return Location(places.cities[0], region, places.countries[0])

# this is faster than geograpy but it still does not work very well
def get_location_geotext(text: str) -> Optional[Location]:
    """extract location using geotext"""
    places = geotext.GeoText(text)
    if not places.countries or not places.cities:
        return None
    return Location(places.cities[0], None, places.countries[0])

def get_papers_with_locations(papers: list[dict[str, Any]]) -> list[Paper]:
    """return list of papers containing list of locations"""
    result: list[Paper] = []
    fail_count = 0
    affiliation_count = 0
    # for paper in progressbar(papers, "extracting locations: "):
    for paper in papers:
        if 'AD' in paper and 'PMID' in paper:
            locations = []
            for affiliation in paper['AD']:
                affiliation_count += 1
                location = get_location_naive(affiliation.lower())
                if location is not None:
                    locations.append(location)
                else:
                    fail_count += 1
            result.append(Paper(int(paper['PMID']), locations))
    if affiliation_count > 0:
        print(f"Extraction fail rate: {round(fail_count / affiliation_count * 100, 2)}% \
(failed: {fail_count}, total: {affiliation_count})")
    return result
