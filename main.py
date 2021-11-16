from plotly.subplots import make_subplots
from processing.searcher import Searcher
import plotly.graph_objects as go
import pandas as pd
import numpy as np


topic = "Alzheimer"
searcher = Searcher(papers="/home/karlo/Desktop/UNI/9. Semester/AI/IntelligenteSysteme/IntelligenteSystemeProject/Collaboration_Network/data/papers/papers.pkl", countryJson="/home/karlo/Desktop/UNI/9. Semester/AI/IntelligenteSysteme/IntelligenteSystemeProject/Collaboration_Network/data/dictionaries/countries.json", topic=topic)
searcher.findTotalResearchPerCountry()
searcher.findTopicInPapers()

countries = []
founds = []
totals = []
ratio = []

for country in searcher.result:
    countries.append(country)
    found = searcher.result[country]["found"]
    total = searcher.result[country]["total"]
    founds.append(found)
    totals.append(total)
    if total > 0:
        ratio.append(round((found/total), ndigits=2))
    else:
        ratio.append(0.0)



costumData = np.empty(shape=(len(totals), 3, 1), dtype='object')
costumData[:, 0] = np.array(totals).reshape(-1, 1)
costumData[:, 1] = np.array(founds).reshape(-1, 1)
costumData[:, 2] = np.array(countries).reshape(-1, 1)

ratioMap = go.Figure(data=go.Choropleth(
    locations=countries,
    z=ratio,
    text=countries,
    locationmode="country names",
    customdata=costumData,
    hovertemplate='<br>Name:%{customdata[2]}<br>ratio:%{z}<br>total:%{customdata[0]}<br>found:%{customdata[1]}',
    colorscale='Blues',
    autocolorscale=False,
    reversescale=False,
    marker_line_color='darkgray',
    marker_line_width=0.5,
    colorbar_title='Ratio of conducted research on the topic ' + topic,
))


ratioMap.update_layout(
    title_text='Feature Karl',
    geo=dict(
        showframe=False,
        showcoastlines=False,
        projection_type='equirectangular'
    ),
    annotations=[dict(
        x=0.55,
        y=0.1,
        xref='paper',
        yref='paper',
        showarrow=False
    )]
)

ratioMap.show()
