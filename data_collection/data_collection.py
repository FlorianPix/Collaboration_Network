import pickle

from Bio import Entrez, Medline


def get_papers(my_query, max_papers, my_email="florian.pix@mailbox.tu-dresden.de"):
    # Get articles from PubMed
    Entrez.email = my_email
    record = Entrez.read(Entrez.esearch(db="pubmed", term=my_query, retmax=max_papers))
    id_list = record["IdList"]
    print("\nThere are %d records for %s."%(len(id_list), my_query.strip()))
    records = Medline.parse(Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text"))
    # records is iterable, which means that it can be consumed only once.
    # Converting it to a list, makes it permanently accessible.
    return list(records)


my_query = "clustering[ti] algorithm"  # query in title and abstract
max_papers = 10  # limit the number of
papers = get_papers(my_query, max_papers)

with open('papers.pkl', 'wb') as outp:
    pickle.dump(papers, outp, pickle.HIGHEST_PROTOCOL)
