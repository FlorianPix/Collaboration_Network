from Bio import Entrez, Medline
from typing import Any
import file_io as io

EFETCH_MAX = 1000


def getPapers(query, max_papers):
    """retrieve papers from PupMed"""
    Entrez.email = "karl_christian.lautenschlaeger@mailbox.tu-dresden.de"
    record = Entrez.read(Entrez.esearch(
        db="pubmed",
        term=query,
        usehistory='y',
        retmax=max_papers if max_papers is not None else 250000))
    id_list = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    print("\nThere are %d papers for %s." % (len(id_list), query.strip()))

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
        result += records
    print(f"number of papers: {len(result)}")
    return result


if __name__ == '__main__':
    paper_dicts = getPapers(query="AD", max_papers=1000)
    io.save_papers_pickle(papers=paper_dicts, filename="papers1000.pkl")
