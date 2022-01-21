from processing.searcher import Searcher
from processing.plotter import Plotter
import plotly.graph_objects as go
import numpy as np

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