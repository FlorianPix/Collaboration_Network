from geotext import GeoText
import time

import pmData, geoData, plotData
    
def main():
    #pmData.getStoreData("pig", 5000)
    
    print("Reading affiliations file... ", end="")
    time0 = time.time()
    affiliationList = pmData.getAllAffiliations()
    time1 = time.time()
    print("done, %s seconds elapsed" %(time1 - time0))
    print("Finding potential cities... ", end="")
    time0 = time.time()
    places = GeoText(" ".join(affiliationList))
    cities = list(places.cities) # list of cities
    time1 = time.time()
    print("done, %s seconds elapsed" %(time1 - time0))
    
    geoData.geolocate(cities)
    plotData.plotCities()

    print("Done")

if __name__ == "__main__":
    main()