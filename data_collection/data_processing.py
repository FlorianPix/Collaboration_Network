import pickle
import json
from typing import Any, Optional
from pprint import pprint

def get_countries() -> dict[str, dict[str, list[str]]]:
    with open('data_collection/countries.json', 'r') as countries_file:
        return json.load(countries_file)

def get_saved_papers(filename) -> list[dict[str, Any]]:
    with open(f'saved_papers/{filename}', 'rb') as papers_file:
        return pickle.load(papers_file)


class Affiliation:
    def __init__(self, front: str, city: str, state: Optional[str], country: str):
        self.front = front
        self.city = city
        self.state = state
        self.country = country


class Paper:
    def __init__(self, id: int, affiliations: list[Affiliation]):
        self.id: int = id
        self.affiliations: list[Affiliation] = affiliations

    def __repr__(self) -> str:
        return f"Paper({self.id}, {len(self.affiliations)} affiliations)"


def extract_country(countries: dict[str, list[str]], affiliation: list[str]):
    target = affiliation[-1].strip()
    for country in countries:
        if country.startswith(target) or target.startswith(country):
            return country
    for country, alt_dict in countries.items():
        for alt in alt_dict["alt"]:
            if alt.startswith(target) or target.startswith(alt):
                return country
    for country, alt_dict in countries.items():
        for code in alt_dict["code"]:
            if target.startswith(code) or target.endswith(code):
                return country
    print(f"Unknown country: '{target}'")
    return None

def get_paper_objects() -> list[Paper]:
    countries = get_countries()
    papers = get_saved_papers('papers.pkl')
    result: list[Paper] = []
    for paper in papers[:100]:
        id = paper.get('PMID', -1)
        affiliations = paper.get('AD', [])
        paper_obj = Paper(id, [])
        for affiliation in affiliations:
            affiliation: str = affiliation
            if affiliation.endswith('.'):
                affiliation = affiliation[:-2]
            parts = affiliation.split(',')
            if len(parts) < 3:
                print("Affiliation too short:", affiliation)
                continue
            country = extract_country(countries, parts)
            if country == None:
                continue
            # TODO: special handling for USA states
            city = parts[-2].strip()
            front = parts[:-2]
            paper_obj.affiliations.append(
                Affiliation(front, city, None, country))
        result.append(paper_obj)
    return result

def main():
    papers = get_paper_objects();
    pprint(papers)

if __name__ == '__main__':
    main()