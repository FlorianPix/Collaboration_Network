import json
import file_io as io
from processing.locator import Locator

"""
Searcher Class:
Counts the number of papers per country and research topic
"""
class Searcher:
    def __init__(self, countryJson, papers, topic: str):
        """
        initialize Searcher class
        coutnryJson: Json, which contains all countries and there coordinates (longitude/ latitude)
        papers: Pickle, which contains a list of python objects that hold all necessary information of a paper
        topic: research topic to which to look for when iterating over the papers
        """
        self.__countries = json.load(open(countryJson, 'rb'))
        self.result = {}
        self.setCountries()
        self.papers = io.get_papers_pickle(papers)
        self.locator = Locator(countries=self.__countries)
        self.topic = ""
        self.setTopic(topic)

    def setCountries(self):
        """
        initialize result dictionary with the following format:
        self.result[country]:
           -total: total number of papers found for that country
           -found: number of papers found for a specific research topic
           -code: country-code
        """
        for country in self.__countries:
            country_result = {"total": 0, "found": 0, "code": self.__countries[country]["code"]}
            self.result[country] = country_result

    def setTopic(self, topic):
        """setter for research topic"""
        self.topic = topic

    def findTotalResearchPerCountry(self):
        """searches for the total number of papers per country"""
        for paper in self.papers:
            if 'AD' in paper and 'PMID' in paper:
                affiliations = paper["AD"]
                foundCountries = self.locator.get_paper_locations(affiliations=affiliations)
                for foundCountry in foundCountries:
                    self.result[foundCountry]["total"] += 1

    def __findTopicInTitle(self, title: str) -> bool:
        """ helper method to decide whether topic is part of title """
        if self.topic.lower() in title.lower():
            return True
        else:
            return False

    def __findTopicInAbstract(self, abstract: str) -> bool:
        """ helper method to decide whether topic is part of abstract """
        if self.topic.lower() in abstract.lower():
            return True
        else:
            return False

    def findTopicInPapers(self):
        """searches for the amount of research in a specific research field per country"""
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


