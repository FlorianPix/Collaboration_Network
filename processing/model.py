"""data classes for model entities"""
from typing import Optional
from dataclasses import dataclass


@dataclass(frozen=True, unsafe_hash=True)
class Location:
    """Class containing a city, a country and optionally a state"""
    city: str
    state: Optional[str]
    country: str
    def __repr__(self):
        return "Location_" + ", ".join([x for x in [self.city, self.state, self.country] if x is not None])
    def __str__(self):
        return ", ".join([x for x in [self.city, self.state, self.country] if x is not None])



@dataclass(frozen=True, unsafe_hash=True)
class Coordinates:
    """Coordinates represented by latitude and longitude (float values)"""
    lat: float
    long: float


@dataclass(frozen=True, unsafe_hash=True)
class Paper:
    """Class representing paper"""
    pmid: int
    locations: list[Location]
