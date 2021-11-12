# Collaboration_Network

## Assignment

Build a network of co-occurring cities in publication affiliations. Visualise the network and find interesting relationships. Which are the strongest co-occurrences across countries? Does geography play a role? Use relative optionally log-odds ratios in your analysis. Draw the network on a map.

## How to use

```bash
python3 -m venv .venv # create virtual environment
source .venv/bin/activate # activate virtual environment
pip install -r requirements.txt # install dependencies
python data_collection.py # to fetch some papers and coordinates
python main.py # to build graph, run visualization, ...
```
* papers and coordinates fetched from api can be saved locally with methods in `file_io.py` (files can be found in `data/papers` and `data/coordinates`)
* a list with countries names, coordinates, codes and alternative names can be found in `data/dictionaries`

## Questions

* How should we use log odds ratios? What is meant by "relative optionally"?
* Concerning "interesting relationships" in the assignment: One idea we had was finding out which countries have a particularly well research in a specific topic. Are we on the right track? Probably that would be a use case for log odds ratios?
* Is it crucial to have great performance when analyzing the data? How big should the dataset be?
* How much accuracy is expected when extracting city names? Probably there will always be affiliations where we won't be able to extract a city and country. (Affiliations don't have a uniform format.)
* (How to fetch **all** papers from PubMed at once?)


# TODO
* [ ] improve city visualization
* [ ] visualize relationships between countries
* [ ] improve location extraction
* [ ] make coords fetching faster (maybe google maps API?)
* [ ] try to find interesting relationships (e. g. per research field)