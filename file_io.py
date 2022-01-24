import pickle

"""
Pickle-file processing:
- saves the downloaded papers into a pickle-file
- load the saved papers from a pickle-file
"""


def save_papers_pickle(papers, filename: str):
    """save papers (dict format) to file in folder data/papers"""
    with open(r"./data/papers/{0}".format(filename), 'wb') as papers_file:
        pickle.dump(papers, papers_file, pickle.HIGHEST_PROTOCOL)


def get_papers_pickle(filePath: str):
    """retrieve papers (dict format) from file in folder data/papers"""
    with open(filePath, 'rb') as papers_file:
        return pickle.load(papers_file)
