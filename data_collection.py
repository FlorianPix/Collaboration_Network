
from io import save_papers_pickle
from typing import Any
from Bio import Entrez, Medline


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


if __name__ == '__main__':
    papers = get_papers("[AD]", None)
    save_papers_pickle(papers, "papers.pkl")
