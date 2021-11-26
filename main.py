from processing.searcher import Searcher
from processing.plotter import Plotter
import plotly.graph_objects as go
import numpy as np


topic = "Alzheimer"
searcher = Searcher(papers="./data/papers/papers.pkl", countryJson="./data/dictionaries/countries.json", topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()

plotter = Plotter(countryPopulation="./data/dictionaries/countryPopulation.json")

plotter.setData(data=searcher.result, topic=topic, paperCount=200000)

plotter.plotResearchRatio()
plotter.plotResearchPopulationRatio()
