
import file_io as io
from processing.extraction import get_papers_with_locations
import processing.network as n
from visualization.interactive_sphere_projection import vis


# run data_collection.py script once to get the files

papers = get_papers_with_locations(io.get_papers_pickle("papers500.pkl"))
coords = io.get_coords_pickle("coordinates500.pkl")

# build a graph with locations (containing a city and country, respectively);
# this graph will have the number of co-occurrences in papers as edge weights
graph = n.build_city_graph(papers, coords)

vis(coords, graph.nodes, graph.edges)
