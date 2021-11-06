import time
from processing.extraction import get_location_naive, get_papers_with_locations
from processing.paper_reading import get_papers_pickle
from processing.coordinates import get_location_coordinates


def test_fail_rate(affiliation_list: list[str]):
    """Print fail rate"""
    fail_counter = 0
    time1 = time.time()
    for affil in affiliation_list:
        # if get_location_geograpy(affiliation) is None:
        # if get_location_geotext(affiliation) is None:
        if get_location_naive(affil) is None:
            fail_counter += 1
    print(time.time() - time1)
    print(f"number of affiliations: {len(affiliation_list)}, fails: {fail_counter}, ratio: {fail_counter / len(affiliation_list)}")


papers = get_papers_pickle("papers010000.pkl")
print(f"number of papers: {len(papers)}")

papers = get_papers_with_locations(papers)
time1 = time.time()
coords = get_location_coordinates(papers)
with open("data/")
print(time.time() - time1)