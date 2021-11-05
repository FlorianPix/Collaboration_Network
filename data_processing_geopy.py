import time
from processing.extraction import get_location_geograpy, get_location_geotext, get_location_naive
from processing.paper_reading import get_papers_pickle

papers = get_papers_pickle("papers010000.pkl")
affiliations = [affiliation for paper in papers for affiliation in paper.get('AD', [])]

fail_counter = 0
time1 = time.time()
for affiliation in affiliations:
    # if get_location_geograpy(affiliation) is None:
    # if get_location_geotext(affiliation) is None:
    if get_location_naive(affiliation) is None:
        fail_counter += 1
print(time.time() - time1)
print(f"number of affiliations: {len(affiliations)}, fails: {fail_counter}, ratio: {fail_counter / len(affiliations)}")
print(f"number of papers: {len(papers)}")
