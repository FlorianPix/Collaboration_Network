
import time
from processing.extraction import get_location_naive

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
    print(f"number of affiliations: {len(affiliation_list)}, "\
        f"fails: {fail_counter}, ratio: {fail_counter / len(affiliation_list)}")
