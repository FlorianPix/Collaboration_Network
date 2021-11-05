from geotext import GeoText
import time

import pmData, geoData, plotData
    
def main():
    #pmData.getStoreData("cancer", 20)
    
    print("Reading affiliations file...")
    time0 = time.time()
    affiliationList = pmData.getAllAffiliations()
    places = GeoText(" ".join(affiliationList))
    time1 = time.time()
    print("Done, %s seconds elapsed" %(time1 - time0))
    cities = list(places.cities) # list of cities
    
    #geoData.geolocate(cities)
    plotData.plotCities()

    print("Done")

if __name__ == "__main__":
    main()