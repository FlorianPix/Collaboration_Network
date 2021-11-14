from Bio import Entrez, Medline


def getPapers(myQuery, maxPapers, myEmail="karl.1998@web.de"):
    # Get articles from PubMed
    Entrez.email = myEmail
    record = Entrez.read(Entrez.esearch(db="pubmed", term=myQuery, retmax=maxPapers))
    idlist = record["IdList"]
    print("\nThere are %d records for %s." % (len(idlist), myQuery.strip()))
    records = Medline.parse(Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text"))
    # records is iterable, which means that it can be consumed only once.
    # Converting it to a list, makes it permanently accessible.
    return list(records)


myQuery = "cancer " + "[tiab]"  # query in title and abstract
maxPapers = 10  # limit the number of papers retrieved
records = getPapers(myQuery, maxPapers)

locationAndPB = []
countOccurrences = []
for record in records:
    plAndPb = (record["PL"], record["PB"])
    if plAndPb not in locationAndPB:
        locationAndPB.append(plAndPb)
        countOccurrences.append(1)
    else:
        index = locationAndPB.index(plAndPb)
        countOccurrences[index] = countOccurrences[index] + 1

#DOI or all AIDs

#not all have references for PB and PL
