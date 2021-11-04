import pickle

from Bio import Entrez, Medline


def get_papers(my_query, max_papers, my_email="oliver.portee@mailbox.tu-dresden.de"):
    # Get articles from PubMed
    Entrez.email = my_email
    record = Entrez.read(Entrez.esearch(db="pubmed", term=my_query, retmax=max_papers))
    id_list = record["IdList"]
    print("\nThere are %d records for %s."%(len(id_list), my_query.strip()))
    records = Medline.parse(Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text"))
    # records is iterable, which means that it can be consumed only once.
    # Converting it to a list, makes it permanently accessible.
    return list(records)

def save_papers(filename, max_number, search_term="[AD]"):
    papers = get_papers(search_term, max_number)
    with open(f'saved_papers/{filename}', 'wb') as outp:
        pickle.dump(papers, outp, pickle.HIGHEST_PROTOCOL)
        print(f'saved {len(papers)} papers in saved_papers/{filename}')

if __name__ == '__main__':
    save_papers('papers.pkl', 100000)