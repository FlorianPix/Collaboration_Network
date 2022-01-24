"""network functions"""
import itertools
from typing import Optional

import networkx as nx

from processing.model import Coordinates, Location, Paper
from processing.util import progressbar


def build_city_graph(papers: list[Paper], coords: dict[Location, Optional[Coordinates]]) \
    -> nx.Graph:
    """
    create graph of locations that have coordinates with number of co-occurrences as edge weights
    """
    graph = nx.Graph()
    for paper in progressbar(papers, "building graph: "):
        valid_locations = filter(lambda l: l in coords and coords[l] is not None, list(set(paper.locations)))
        for (l_1, l_2) in itertools.combinations(valid_locations, 2):
            if graph.has_edge(l_1, l_2):
                graph[l_1][l_2]['weight'] += 1
            else:
                graph.add_edge(l_1, l_2, weight=1)
    return graph

def build_country_graph(papers: list[Paper], coords: dict[Location, Optional[Coordinates]]) -> nx.Graph:
    """
    create graph of locations with number of co-occurrences as edge weights
    """
    graph = nx.Graph()
    for paper in progressbar(papers, "building graph: "):
        valid_locations = filter(lambda l: l in coords and coords[l] is not None, list(set(paper.locations)))
        valid_locations = {Location(city=None, state=None, country=l.country) for l in valid_locations}
        for (l_1, l_2) in itertools.combinations_with_replacement(valid_locations, 2):
            l_1: Location
            l_2: Location
            if graph.has_edge(l_1, l_2):
                graph[l_1][l_2]['weight'] += 1
            else:
                graph.add_edge(l_1, l_2, weight=1)
    return graph
