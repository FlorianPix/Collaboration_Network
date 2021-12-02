"""unpickle papers and coords, build and visualize graph, analyze data"""
import file_io as io

from processing.extraction import get_papers_with_locations
from processing import network as n, analysis as a
from visualization.interactive_sphere_projection import vis

# run data_collection.py script once to get the files
DATASET_SIZE = 10000
papers = get_papers_with_locations(io.get_papers_pickle(f"papers{DATASET_SIZE}.pkl"))
coords = io.get_coords_pickle(f"coordinates{DATASET_SIZE}.pkl")

# build a graph with locations (containing a city and country, respectively);
# this graph will have the number of co-occurrences in papers as edge weights
city_graph = n.build_city_graph(papers, coords)
country_graph = n.build_country_graph(papers, coords)

# filter countries by paper threshold and sort alphabetically
THRESHOLD = 10
country_pubs = a.publications_by_country(papers)
relevant_countries = list(filter(lambda x:country_pubs[x]>=THRESHOLD, country_graph.nodes))
relevant_countries = list(sorted(relevant_countries, key=lambda x:x.country))

log_odds = a.log_odds_countries(country_graph, papers, country_graph.nodes, THRESHOLD)
out, errors = "", ""
for country in relevant_countries:
    try:
        out += f"Log odds for {country} ({country_pubs[country]} publications):" \
            + f"\t\t{sorted(log_odds[country].items(), key=lambda x:x[1], reverse=True)}\n"
    except Exception as e:
        errors += f"Couldn't determine log odds ratios for {country}: {e!r}\n"
with open(f"output{DATASET_SIZE}_{THRESHOLD}.txt", mode="wt", encoding="utf-8") as output_file:
    output_file.write(out)
with open(f"errors{DATASET_SIZE}_{THRESHOLD}.txt", mode="wt", encoding="utf-8") as errors_file:
    errors_file.write(errors)
