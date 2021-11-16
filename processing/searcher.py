import json
import file_io as io
from processing.locator import Locator


class Searcher:
    def __init__(self, countryJson, papers, topic: str):
        self.__countries = json.load(open(countryJson, 'rb'))
        self.result = {}
        self.setCountries()
        self.papers = io.get_papers_pickle(papers)
        self.locator = Locator(countries=self.__countries)
        self.topic = ""
        self.setTopic(topic)

    def setCountries(self):
        for country in self.__countries:
            country_result = {"total": 0, "found": 0, "code": self.__countries[country]["code"]}
            self.result[country] = country_result

    def setTopic(self, topic):
        self.topic = topic

    def findTotalResearchPerCountry(self):
        for paper in self.papers:
            if 'AD' in paper and 'PMID' in paper:
                affiliations = paper["AD"]
                foundCountries = self.locator.get_paper_locations(affiliations=affiliations)
                for foundCountry in foundCountries:
                    self.result[foundCountry]["total"] += 1

    def __findTopicInTitle(self, title: str) -> bool:
        if self.topic.lower() in title.lower():
            return True
        else:
            return False

    def __findTopicInAbstract(self, abstract: str) -> bool:
        if self.topic.lower() in abstract.lower():
            return True
        else:
            return False

    def findTopicInPapers(self):
        for paper in self.papers:
            if 'AD' in paper and 'PMID' in paper:
                found = False
                if 'TI' in paper:
                    found = self.__findTopicInTitle(title=paper['TI'])
                if not found and 'AB' in paper:
                    found = self.__findTopicInAbstract(abstract=paper['AB'])
                if found:
                    affiliations = paper["AD"]
                    foundCountries = self.locator.get_paper_locations(affiliations=affiliations)
                    for foundCountry in foundCountries:
                        self.result[foundCountry]["found"] += 1


