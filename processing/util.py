"""helper functions"""
import sys

def progressbar(items, prefix="", size=60, file=sys.stdout):
    """display progressbar; source: https://stackoverflow.com/a/34482761"""
    count = len(items)
    def show(j):
        progess = int(size*j/count)
        file.write(f"{prefix}[{'#' * progess}{'.' * (size - progess)}] {j}/{count}\r")
        file.flush()
    show(0)
    for i, item in enumerate(items):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()
