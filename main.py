from processing.searcher import Searcher
from processing.plotter import Plotter
import plotly.graph_objects as go
import numpy as np

paperCount = 200000
topic = "Alzheimer"
searcher = Searcher(papers="/home/karlo/Desktop/UNI/9. Semester/AI/IntelligenteSysteme/IntelligenteSystemeProject/Collaboration_Network/data/papers/papers.pkl", countryJson="/home/karlo/Desktop/UNI/9. Semester/AI/IntelligenteSysteme/IntelligenteSystemeProject/Collaboration_Network/data/dictionaries/countries.json", topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()
plotter = Plotter(countryPopulation="/home/karlo/Desktop/UNI/9. Semester/AI/IntelligenteSysteme/IntelligenteSystemeProject/Collaboration_Network/data/dictionaries/countryPopulation.json")
plotter.setData(data=searcher.result, topic=topic, paperCount=paperCount)
plotter.plotResearchRatio()
plotter.plotResearchPopulationRatio()

# new topic:
topic = "Tourette"
searcher.setCountries()
searcher.setTopic(topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()
plotter.setData(data=searcher.result, topic=topic, paperCount=paperCount)
plotter.plotResearchRatio()