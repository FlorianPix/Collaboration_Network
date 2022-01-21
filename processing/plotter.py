import json
import plotly.graph_objects as go
import numpy as np


class Plotter:
    def __init__(self, countryPopulation):
        self.__countryPopulation = json.load(open(countryPopulation, 'rb'))
        self.countries = []
        self.founds = []
        self.totals = []
        self.ratioResearch = []
        self.ratioPopulation = []
        self.population = []
        self.topic = ""
        self.paperCount = 0

    def __resetData(self):
        self.countries = []
        self.founds = []
        self.totals = []
        self.ratioResearch = []
        self.ratioPopulation = []
        self.population = []
        self.topic = ""
        self.paperCount = 0

    def setData(self, data, topic, paperCount):
        self.__resetData()
        for country in data:
            self.countries.append(country)
            found = data[country]["found"]
            total = data[country]["total"]
            self.founds.append(found)
            self.totals.append(total)


            if country in self.__countryPopulation:
                population = self.__countryPopulation[country]["population"]
                self.population.append(population)
            else:
                self.population.append(0)
            if total > 0:
                self.ratioResearch.append(round((found / total), ndigits=2))
                if population > 0:
                    self.ratioPopulation.append(round((total / population), ndigits=10))
                else:
                    self.ratioPopulation.append(0.0)
            else:
                self.ratioResearch.append(0.0)
                self.ratioPopulation.append(0.0)

        # convert list to numpy array for better handling
        ratioPopNp = np.array(self.ratioPopulation)

        # get the minimum value of the array which is non-zero
        minPopulationRatio = np.min(np.where(ratioPopNp>0,ratioPopNp,np.max(ratioPopNp)))

        # adjust the value space of the ratios to be y>=1 for positive values in log-space
        ratioPopNp = np.where(ratioPopNp>0,ratioPopNp/minPopulationRatio,1)

        # calculate the values in log10 space (unit: dB)
        ratioPopLog = 10*np.log10(ratioPopNp)

        # get original data format (list)
        self.ratioPopulation = ratioPopLog.tolist()


        self.topic = topic
        self.paperCount = paperCount

    def plotResearchRatio(self):
        customData = np.empty(shape=(len(self.totals), 3, 1), dtype='object')
        customData[:, 0] = np.array(self.totals).reshape(-1, 1)
        customData[:, 1] = np.array(self.founds).reshape(-1, 1)
        customData[:, 2] = np.array(self.countries).reshape(-1, 1)
        researchRatioMap = go.Figure(data=go.Choropleth(
            locations=self.countries,
            z=self.ratioResearch,
            text=self.countries,
            locationmode="country names",
            customdata=customData,
            hovertemplate='<br>Name:%{customdata[2]}<br>ratio:%{z}<br>total:%{customdata[0]}<br>found:%{customdata[1]}',
            colorscale='Blues',
            autocolorscale=False,
            reversescale=False,
            marker_line_color='darkgray',
            marker_line_width=0.5,
            colorbar_title='Ratio of medical research on the topic {0} to a country\'s total medical research'.format(
                self.topic),
        ))
        researchRatioMap.update_layout(
            title_text='Map showing which countries are focusing on the research topic {0}'.format(self.topic),
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
                showarrow=False,
                text="Number of analysed papers: {0}".format(self.paperCount)
            )]
        )

        researchRatioMap.show()

    def plotResearchPopulationRatio(self):
        customData = np.empty(shape=(len(self.totals), 3, 1), dtype='object')
        customData[:, 0] = np.array(self.totals).reshape(-1, 1)
        customData[:, 1] = np.array(self.population).reshape(-1, 1)
        customData[:, 2] = np.array(self.countries).reshape(-1, 1)
        researchRatioMap = go.Figure(data=go.Choropleth(
            locations=self.countries,
            z=self.ratioPopulation,
            text=self.countries,
            locationmode="country names",
            customdata=customData,
            hovertemplate='<br>Name:%{customdata[2]}<br>ratio:%{z}<br>total:%{customdata[0]}<br>population:%{customdata[1]}',
            colorscale='Blues',
            autocolorscale=False,
            reversescale=False,
            marker_line_color='darkgray',
            marker_line_width=0.5,
            colorbar_title='Ratio of public medical research papers to a country\'s population [dB]',
        ))
        researchRatioMap.update_layout(
            title_text='Map, showing which countries do a lot of medical research wrt to population of the country',
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
                showarrow=False,
                text="Number of analysed papers: {0}".format(self.paperCount)
            )]
        )

        researchRatioMap.show()


