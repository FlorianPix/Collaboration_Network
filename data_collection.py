import pickle
from typing import Any
from Bio import Entrez, Medline
from pprint import pprint


EFETCH_MAX = 10000


def get_papers(my_query, max_papers) -> list[dict[str, Any]]:
    # Get articles from PubMed
    Entrez.email = "oliver.portee@mailbox.tu-dresden.de"
    record = Entrez.read(Entrez.esearch(
        db="pubmed",
        term=my_query,
        usehistory='y',
        retmax=max_papers if max_papers is not None else 250000))
    id_list = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    print("\nThere are %d records for %s." % (len(id_list), my_query.strip()))

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


def save_papers(filename, max_number, search_term="[AD]"):
    papers = get_papers(search_term, max_number)
    # pprint(list(map(lambda p: p['PMID'], papers)))
    with open(f'papers/{filename}', 'wb') as outp:
        pickle.dump(papers, outp, pickle.HIGHEST_PROTOCOL)
        print(f'saved {len(papers)} papers in data/papers/{filename}')


if __name__ == '__main__':
    save_papers('papers.pkl', None)
