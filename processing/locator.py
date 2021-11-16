from typing import Any, Optional


class Locator:
    def __init__(self, countries):
        self.__countries = countries

    def __extract_country(self, affiliation) -> Optional[str]:
        target = affiliation[-1].strip()
        for country in self.__countries:
            if country.startswith(target) or target.startswith(country):
                return country
        for country, data in self.__countries.items():
            for alt in data["alt"]:
                if alt.startswith(target) or target.startswith(alt):
                    return country
        for country, data in self.__countries.items():
            code = data["code"]
            if target.startswith(code) or target.endswith(code):
                return country
        # print(f"Unknown country: '{target}'")
        return None

    # there's a lot of room for optimization here
    def get_location_naive(self, text: str) -> Optional[str]:
        """extract location using naive parsing"""
        if text.endswith('.'):
            text = text[:-2]
        parts = text.split(',')
        if len(parts) < 3:
            # print("Affiliation too short:", text)
            return None
        country = self.__extract_country(parts)
        if country is None:
            return None
        # TODO: special handling for USA states
        return country

    def get_paper_locations(self, affiliations):
        """return list of papers containing list of locations"""
        foundCountries = []
        for affiliation in affiliations:
            country = self.get_location_naive(text=affiliation)
            if country is not None and country not in foundCountries:
                foundCountries.append(country)
        return foundCountries
