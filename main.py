from processing.searcher import Searcher
from processing.plotter import Plotter

"""
Main:
- creates objects of searcher and plotter and executes their functionalities
- the script searches for papers on "Alzheimer" and "cancer" and plots the ratio of found papers on the specified topic
 in regards to the total amount of research per country
- the script also searches for the total amount of papers published in an country and plots the ratio of the
 found papers in regards to the total population of the country.

"""

paperCount = 250000
topic = "Alzheimer"
searcher = Searcher(papers="./data/papers/papers250000.pkl", countryJson="./data/dictionaries/countries.json", topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()
plotter = Plotter(countryPopulation="./data/dictionaries/countryPopulation.json")
plotter.setData(data=searcher.result, topic=topic, paperCount=paperCount)
plotter.plotResearchRatio()
plotter.plotResearchPopulationRatio()

# new topic:
topic = "cancer"
searcher.setCountries()
searcher.setTopic(topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()
plotter.setData(data=searcher.result, topic=topic, paperCount=paperCount)
plotter.plotResearchRatio()