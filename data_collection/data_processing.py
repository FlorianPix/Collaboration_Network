import pickle
import wikipedia
import json

f = open("countries.json")
countries = json.loads(f.read())

"""
links = wikipedia.page("list of cities by country").links
links_to_lists_of_cities_by_country = list()
for link in links:
    if link.startswith("List of"):
        links_to_lists_of_cities_by_country.append(link)
        print(link)
"""

papers = pickle.load(open('papers.pkl', 'rb'))

for p in papers:
    print(p['PL'])
    try:
        ad = p['AD']
        for a in ad:
            sp = a.split(',')
            for s in sp:
                s = s.replace(".", "")
    except:
        print('no department')
