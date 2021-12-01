"""unpickle papers and coords, build and visualize graph"""

import json
import file_io as io

from processing.extraction import get_papers_with_locations
from processing import network as n, analysis as a
from visualization.interactive_sphere_projection import vis

# run data_collection.py script once to get the files

papers = get_papers_with_locations(io.get_papers_pickle("papers100.pkl"))
coords = io.get_coords_pickle("coordinates100.pkl")

# build a graph with locations (containing a city and country, respectively);
# this graph will have the number of co-occurrences in papers as edge weights
graph = n.build_city_graph(papers, coords)

vis(coords, graph.nodes, graph.edges(data=True))

#print(a.publications_per_country(papers))
dic = a.country_collaborations(graph)
for c in json.load(open("data/dictionaries/countries.json", "rb")):
    input(c + str({k: v for k, v in (sorted(dic[c].items(), key=lambda x:x[1], reverse=True)) if v > 0}))
