from Bio import Entrez, Medline
import json, time

def downloadPapers(query, maxPapers, email = "johannalbrecht@protonmail.com"):
    Entrez.email = email
    print("Accessing PM database...")
    time0 = time.time()
    record = Entrez.read(Entrez.esearch(db="pubmed", term=query+"[tiab]", retmax=maxPapers))
    idlist =record["IdList"]
    print("There are %d records for %s." %(len(idlist), query.strip()))
    records = Medline.parse(Entrez.efetch(db="pubmed",id=idlist, rettype="medline", \
        retmode="text"))
    time1 = time.time()
    print("Done, %s seconds elapsed" %(time1 - time0))
    # records is iterable, which means that it can be consumed only once.
    # Converting it to a list, makes it permanently accessible.
    return list(records)

def getStoreData(query, maxPapers):
    records = downloadPapers(query, maxPapers)
    affils = []
    print("Writing to file...")
    time0 = time.time()
    for record in records:
        if record.get("AD"):
            affils.append(" ".join(record.get("AD")))
    with open("affiliations.json", "w") as file:
        json.dump(affils, file, indent=4)
    time1 = time.time()
    print("Done, %s seconds elapsed" %(time1 - time0))

def getAllAffiliations():
    with open("affiliations.json", "r") as file:
        return json.load(file)