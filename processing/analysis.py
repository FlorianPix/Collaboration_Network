"""testing some analysis of the data"""
import itertools
import json
import math
from typing import Any

import networkx as nx

from processing.model import Location

def publications_by_country(papers: dict[str, Any]) -> dict[Location, int]:
    """returns number of published papers per country"""
    countries_publications = {}
    for country in json.load(open("data/dictionaries/countries.json", "rb")):
        countries_publications[Location(city=None, state=None, country=country)] = 0
    for paper in papers:
        participant_countries = {Location(city=None, state=None, country=location.country) \
            for location in paper.locations}
        for country in participant_countries:
            try:
                countries_publications[country] += 1
            except KeyError:
                countries_publications[country] = 1
    return (dict(sorted(countries_publications.items(), key=lambda x: x[1], reverse=True)))

def log_odds_countries(graph: nx.Graph, papers: dict[str, Any], countries: list[Location], \
    threshold=0) -> nx.Graph:
    """log odds of two countries collaborating"""
    pub_numbers = publications_by_country(papers)
    log_odds_graph = nx.Graph()
    # ignore countries that published less papers than the threshold value
    countries = list(filter(lambda x: pub_numbers[x]>=threshold, countries))
    pubs_total = len(papers)
    for country1, country2 in itertools.combinations_with_replacement(countries, 2):
        if graph.has_edge(country1, country2):
            pubs_overlap = graph[country1][country2]["weight"]
            pubs_country1 = pub_numbers[country1]
            pubs_country2 = pub_numbers[country2]
            log_odds = math.log2(pubs_total*pubs_overlap / (pubs_country1*pubs_country2))
        else:
            log_odds = -math.inf
        log_odds_graph.add_edge(country1, country2, weight=log_odds)
    return log_odds_graph
