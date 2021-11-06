
import pickle
from typing import Any, Optional

from processing.model import Coordinates, Location

def save_papers_pickle(papers: list[dict[str, Any]], filename: str):
    """save papers (dict format) to file in folder data/papers"""
    with open(f"data/papers/{filename}", 'wb') as papers_file:
        pickle.dump(papers, papers_file, pickle.HIGHEST_PROTOCOL)

def get_papers_pickle(filename: str) -> dict[str, Any]:
    """retrive papers (dict format) from file in folder data/papers"""
    with open(f"data/papers/{filename}", 'rb') as papers_file:
        return pickle.load(papers_file)

def save_coords_pickle(coords: dict[Location, Optional[Coordinates]], filename: str):
    """save coordinates (dict with locations as keys) to file in folder data/coordinates"""
    with open(f"data/coordinates/{filename}", 'wb') as coords_file:
        pickle.dump(coords, coords_file, pickle.HIGHEST_PROTOCOL)

def get_coords_pickle(filename: str) -> dict[Location, Optional[Coordinates]]:
    """retrieve coordinates (dict with locations as keys) from file in folder data/coordinates"""
    with open(f"data/coordinates/{filename}", 'rb') as coords_file:
        return pickle.load(coords_file)
