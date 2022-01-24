"""unpickle papers and coords, build and visualize graph, analyze data"""

import json
import math

import file_io as io
import processing.network as n
from processing import analysis as a
from processing.extraction import get_papers_with_locations
from processing.model import Coordinates, Location
from visualization.interactive_sphere_projection import vis, vis2

# run data_collection.py script once to get the files

DATASET_SIZE = 100000
papers = get_papers_with_locations(io.get_papers_pickle(f"papers{DATASET_SIZE}.pkl"))
coords = io.get_coords_pickle(f"coordinates{DATASET_SIZE}.pkl")
countries = json.load(open("data/dictionaries/countries.json", "rb"))
for country in countries.items():
    coords[Location(city=None, state=None, country=country[0])] \
        = Coordinates(lat=country[1]["lat"], long=country[1]["long"])

# build a graph with locations (containing a city and country, respectively);
# this graph will have the number of co-occurrences in papers as edge weights
graph = n.build_city_graph(papers, coords)
country_graph = n.build_country_graph(papers, coords)

# vis(coords, graph.nodes, graph.edges) # city wise
vis2(graph.nodes, graph.edges) # country wise

# minimum number of papers a city/country must have published to be incorporated
THRESHOLD = 20
log_odds_countrygraph = a.log_odds_countries(country_graph, papers, country_graph.nodes, THRESHOLD)
log_odds_citygraph = a.log_odds_cities(graph, papers, graph.nodes, THRESHOLD)

# print into readable format
out = ""
for key, value in a.publications_by_country(papers).items():
    out += f"{key}: {value}\n"
with open(f"publications{DATASET_SIZE}.txt", mode="wt", encoding="utf-8") as output_file:
    output_file.write(out)

out = ""
for key, value in a.publications_by_city(papers).items():
    if value >= 2:
        out += f"{key}: {value}\n"
with open(f"publicationscities{DATASET_SIZE}.txt", mode="wt", encoding="utf-8") as output_file:
    output_file.write(out)

out = ""
country_pubs = a.publications_by_country(papers)
for node in list(sorted(log_odds_countrygraph.nodes, key=lambda x:x.country)):
    out += f"Country: {node}\nTotal publications: {country_pubs[node]}\nLog odds ratios:\n"
    for connection in dict(sorted(log_odds_countrygraph[node].items(), \
        key=lambda x:x[1]["weight"], reverse=True)):
        weight = log_odds_countrygraph[node][connection]['weight']
        if weight > -math.inf:
            out += f"\t{connection}: {log_odds_countrygraph[node][connection]['weight']}\n"
    out += "\n"
with open(f"output{DATASET_SIZE}_{THRESHOLD}.txt", mode="wt", encoding="utf-8") as output_file:
    output_file.write(out)

out = ""
city_pubs = a.publications_by_city(papers)
for node in list(sorted(log_odds_citygraph.nodes, key=lambda x:x.city)):
    out += f"City: {node}\nTotal publications: {city_pubs[node]}\nLog odds ratios:\n"
    for connection in dict(sorted(log_odds_citygraph[node].items(), \
        key=lambda x:x[1]["weight"], reverse=True)):
        weight = log_odds_citygraph[node][connection]['weight']
        if weight > -math.inf:
            out += f"\t{connection}: {log_odds_citygraph[node][connection]['weight']}\n"
    out += "\n"
with open(f"outputcities{DATASET_SIZE}_{THRESHOLD}.txt", mode="wt", encoding="utf-8") \
    as output_file:
    output_file.write(out)
