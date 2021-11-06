
import file_io as io
from processing.extraction import get_papers_with_locations
import processing.network as n

papers = get_papers_with_locations(io.get_papers_pickle("papers010000.pkl"))
coords = io.get_coords_pickle("coordinates010000.pkl")

graph = n.build_city_graph(papers, coords)
for location in graph.nodes:
    print(location, coords[location])
print(graph)