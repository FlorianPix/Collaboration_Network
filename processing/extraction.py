"""Extract location data out of string"""
from typing import Any, Optional
import json
import geograpy
import geotext

from .model import Location, Paper

__countries = json.load(open("data/dictionaries/countries.json", 'rb'))


def __extract_country(affiliation: list[str]) -> Optional[str]:
    target = affiliation[-1].strip()
    for country in __countries:
        if country.startswith(target) or target.startswith(country):
            return country
    for country, alt_dict in __countries.items():
        for alt in alt_dict["alt"]:
            if alt.startswith(target) or target.startswith(alt):
                return country
    for country, alt_dict in __countries.items():
        for code in alt_dict["code"]:
            if target.startswith(code) or target.endswith(code):
                return country
    ## print(f"Unknown country: '{target}'")
    return None

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

def get_location_geograpy(text: str) -> Optional[Location]:
    """extract location using geograpy"""
    places = geograpy.get_geoPlace_context(text=text)
    if not places.cities or not places.countries:
        return None
    region = places.regions[0] if places.regions else None
    return Location(places.cities[0], region, places.countries[0])

def get_location_geotext(text: str) -> Optional[Location]:
    """extract location using geotext"""
    places = geotext.GeoText(text)
    if not places.countries or not places.cities:
        return None
    return Location(places.cities[0], None, places.countries[0])

def get_papers_with_locations(papers: list[dict[str, Any]]) -> list[Paper]:
    """return list of papers containing list of locations"""
    result: list[Paper] = []
    for paper in papers:
        if 'AD' in paper and 'PMID' in paper:
            locations = []
            for affiliation in paper['AD']:
                location = get_location_naive(affiliation)
                if location is not None:
                    locations.append(location)
            result.append(Paper(int(paper['PMID']), locations))
    return result
