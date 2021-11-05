
import pickle
from typing import Any

def get_papers_pickle(filename: str) -> dict[str, Any]:
    with open(f"data/saved_papers/{filename}", 'rb') as papers_file:
        return pickle.load(papers_file)