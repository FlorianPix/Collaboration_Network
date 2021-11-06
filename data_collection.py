
from typing import Any
from Bio import Entrez, Medline
import file_io as io
from processing.extraction import get_papers_with_locations
from processing.coordinates import calculate_paper_coordinates


EFETCH_MAX = 10000


def get_papers(query, max_papers) -> list[dict[str, Any]]:
    """retrieve papers from PupMed"""
    Entrez.email = "oliver.portee@mailbox.tu-dresden.de"
    record = Entrez.read(Entrez.esearch(
        db="pubmed",
        term=query,
        usehistory='y',
        retmax=max_papers if max_papers is not None else 250000))
    id_list = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    print("\nThere are %d records for %s." % (len(id_list), query.strip()))

    result: list[dict[str, Any]] = []
    for i in range(0, len(id_list), EFETCH_MAX):
        records = list(Medline.parse(Entrez.efetch(
            db="pubmed",
            id=id_list,
            retstart=i,
            retmax=EFETCH_MAX,
            rettype="medline",
            retmode="text",
            webenv=webenv,
            query_key=query_key)))
        print(f"i: {i}, length: {len(records)}")
        result += records
    print(f"number of papers: {len(result)}")
    return result

def data_preparation(number_of_papers):
    """fetch papers and calculate coordinates; this can take very long :)"""
    suffix = str(number_of_papers)
    # get paper dictionaries from PubMed
    paper_dicts = get_papers("[AD]", number_of_papers)
    # save papers in `data/papers`
    io.save_papers_pickle(paper_dicts, f"papers{suffix}.pkl")
    # get actual paper objects with locations (containing city, maybe a state and country)
    papers_objs = get_papers_with_locations(paper_dicts)
    print("starting to fetch location coordinates")
    # fetch location coordinates from Nominatim API (this will take some time)
    coordinates = calculate_paper_coordinates(papers_objs)
    # save coordinates in `data/coordinates`
    io.save_coords_pickle(coordinates, f"coordinates{suffix}.pkl")

if __name__ == '__main__':
    data_preparation(10000) # fetch 10000 papers
