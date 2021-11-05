from Bio import Entrez, Medline
import json, time

def downloadPapers(query, maxPapers, email = "johannalbrecht@protonmail.com"):
    Entrez.email = email
    print("Accessing PM database...", end=" ")
    time0 = time.time()
    record = Entrez.read(Entrez.esearch(db="pubmed", term=query+"[tiab]", retmax=maxPapers))
    idlist =record["IdList"]
    print("getting %d records for %s..." %(len(idlist), query.strip()), end=" ")
    records = Medline.parse(Entrez.efetch(db="pubmed", id=idlist, rettype="medline", \
        retmode="text"))
    # records is iterable, which means that it can be consumed only once.
    # Converting it to a list, makes it permanently accessible.
    recordsList = list(records)
    time1 = time.time()
    print("done, %s seconds elapsed" %(time1 - time0))
    return recordsList

def getStoreData(query, maxPapers):
    records = downloadPapers(query, maxPapers)
    affils = []
    print("Writing to file...", end=" ")
    time0 = time.time()
    for record in records:
        if record.get("AD"):
            affils.append(" ".join(record.get("AD")))
    with open("data/affiliations.json", "w") as file:
        json.dump(affils, file, indent=4)
    time1 = time.time()
    print("done, %s seconds elapsed" %(time1 - time0))

def getAllAffiliations():
    with open("data/affiliations.json", "r") as file:
        return json.load(file)