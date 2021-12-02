"""testing some analysis of the data"""
import itertools
import json
import math
from typing import Any, Tuple

import networkx
import numpy as np

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

def country_collaborations(graph: networkx.Graph) -> dict[str, dict[str, int]]:
    """no. of publications where two specific countries collaborated"""
    collab_matrix = {}
    # initialize matrix with zeroes
    for country1 in json.load(open("data/dictionaries/countries.json", "rb")):
        collab_matrix[Location(city=None, state=None, country=country1)] = {}
        for country2 in json.load(open("data/dictionaries/countries.json", "rb")):
            collab_matrix[Location(city=None, state=None, country=country1)] \
                [Location(city=None, state=None, country=country2)] = 0
    for edge in graph.edges(data=True):
        edge: Tuple(Location, Location, dict)
        collab_matrix[edge[0]][edge[1]] = edge[2]["weight"]
        collab_matrix[edge[1]][edge[0]] = edge[2]["weight"]
    return collab_matrix

def log_odds_countries(graph: networkx.Graph, papers: dict[str, Any], countries: list[Location], \
    threshold=0) -> dict[Location, dict[Location, int]]:
    """log odds of two countries collaborating"""
    pub_numbers = publications_by_country(papers)
    collabs = country_collaborations(graph)
    log_odds = {}
    # ignore countries that published less papers than the threshold value
    countries = list(filter(lambda x: pub_numbers[x]>=threshold, countries))
    for country1 in countries:
        log_odds[country1] = {}
        for country2 in countries:
            log_odds[country2] = {}
    pubs_total = len(papers)
    for country1, country2 in itertools.combinations_with_replacement(countries, 2):
        pubs_country1 = pub_numbers[country1]
        pubs_country2 = pub_numbers[country2]
        pubs_overlap = collabs[country1][country2]
        try:
            odds = math.log2(pubs_total*pubs_overlap / (pubs_country1*pubs_country2))
        except (ValueError, ZeroDivisionError):
            odds = -math.inf
        log_odds[country1][country2] = odds
        log_odds[country2][country1] = odds
    return log_odds
