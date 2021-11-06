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

* Is it crucial to have great performance? How big should the dataset be?
* How should be use log odds ratios? What is meant by "relative optionally"?
* How much accuracy is expected when extracting city names? (Affiliations don't have uniform format.)
* (How to fetch **all** papers from PubMed at once?)


# TODO
* [ ] visualization on map
* [ ] visualization on globe
* [ ] improve location extraction
* [ ] make coords fetching faster (maybe google maps API?)
* [ ] try to find interesting relationships (e. g. per field)