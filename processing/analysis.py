"""testing some analysis of the data"""
import json
from typing import Any, Tuple

import networkx

from processing.model import Location

def publications_by_country(papers: dict[str, Any]) -> dict[str, int]:
    """returns number of published papers per country"""
    countries_publications = {}
    for p in papers:
        publisher_countries = set()
        for l in p.locations:
            publisher_countries.add(l.country)
        for c in publisher_countries:
            try:
                countries_publications[c] += 1
            except KeyError:
                countries_publications[c] = 1
    return (dict(sorted(countries_publications.items(), key=lambda x: x[1], reverse=True)))

def country_collaborations(graph: networkx.Graph) -> dict[dict[str, int]]:
    """no. of publications where two specific countries collaborated"""
    collab_matrix = {}
    for c in json.load(open("data/dictionaries/countries.json", "rb")):
        collab_matrix[c] = {}
        for c2 in json.load(open("data/dictionaries/countries.json", "rb")):
            collab_matrix[c][c2] = 0
    for e in graph.edges:
        e: Tuple(Location, Location)
        c1 = e[0].country
        c2 = e[1].country
        collab_matrix[c1][c2] += 1
        collab_matrix[c2][c1] += 1
    return collab_matrix
