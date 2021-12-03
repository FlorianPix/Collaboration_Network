"""unpickle papers and coords, build and visualize graph, analyze data"""

import json

import file_io as io

from processing.extraction import get_papers_with_locations
from processing import network as n, analysis as a
from processing.model import Coordinates, Location
from visualization.interactive_sphere_projection import vis

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
city_graph = n.build_city_graph(papers, coords)
country_graph = n.build_country_graph(papers, coords)

# minimum number of papers a country must have published to be incorporated
THRESHOLD = 20
log_odds_graph = a.log_odds_countries(country_graph, papers, country_graph.nodes, THRESHOLD)

# print into readable format
# code needs improvement
country_pubs = a.publications_by_country(papers)
out = ""
for node in list(sorted(log_odds_graph.nodes, key=lambda x:x.country)):
    out += f"Country: {node}\nTotal publications: {country_pubs[node]}\nLog odds ratios:\n"
    for connection in dict(sorted(log_odds_graph[node].items(), \
        key=lambda x:x[1]["weight"], reverse=True)):
        out += f"\t{connection}: {log_odds_graph[node][connection]['weight']}\n"
    out += "\n"
with open(f"output{DATASET_SIZE}_{THRESHOLD}.txt", mode="wt", encoding="utf-8") as output_file:
    output_file.write(out)

vis(coords, log_odds_graph.nodes, list(filter(lambda x:x[2]["weight"] >= 0, \
    log_odds_graph.edges(data=True))), \
        title_text=f"Log odds ratios of countries, threshold: {THRESHOLD}")
