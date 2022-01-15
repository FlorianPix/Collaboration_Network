"""Extract location data out of string"""
from typing import Any, Optional
import time
import json
import geograpy
import geotext

from .model import Location, Paper
from .util import progressbar

__countries = json.load(open("data/dictionaries/countries.json", 'rb'))

def __extract_country(affiliation: list[str]) -> Optional[str]:
    target = affiliation[-1].strip()
    for country in __countries:
        if country.startswith(target) or target.startswith(country):
            return country
    for country, data in __countries.items():
        for alt in data["alt"]:
            if alt.startswith(target) or target.startswith(alt):
                return country
    for country, data in __countries.items():
        code = data["code"]
        if target.startswith(code) or target.endswith(code):
            return country
    ## print(f"Unknown country: '{target}'")
    return None

# there's a lot of room for optimization here
def get_location_naive(text: str) -> Optional[Location]:
    """extract location using naive parsing"""
    if text.endswith('.'):
        text = text[:-2]
    parts = text.split(',')
    if len(parts) < 3:
        ## print("Affiliation too short:", text)
        return None
    country = __extract_country(parts)
    if country is None:
        return None
    # TODO: special handling for USA states
    city = parts[-2].strip()
    return Location(city, None, country)

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
    for paper in progressbar(papers, "extracting locations: "):
        if 'AD' in paper and 'PMID' in paper:
            locations = []
            for affiliation in paper['AD']:
                affiliation_count += 1
                location = get_location_naive(affiliation)
                if location is not None:
                    locations.append(location)
                else:
                    fail_count += 1
            result.append(Paper(int(paper['PMID']), locations))
    if affiliation_count > 0:
        print(f"Extraction fail rate: {round(fail_count / affiliation_count * 100, 2)}% (failed: {fail_count}, total: {affiliation_count})")
    return result

def test_fail_rate(affiliation_list: list[str]):
    """Print fail rate"""
    fail_counter = 0
    time1 = time.time()
    for affil in affiliation_list:
        # if get_location_geograpy(affiliation) is None:
        # if get_location_geotext(affiliation) is None:
        if get_location_naive(affil) is None:
            fail_counter += 1
    print(time.time() - time1)
    print(f"number of affiliations: {len(affiliation_list)}, "\
        f"fails: {fail_counter}, ratio: {fail_counter / len(affiliation_list)}")
