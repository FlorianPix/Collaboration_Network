
from typing import Optional


class Location:
    """Class containing a city, a country and optionally a state"""

    def __init__(self, city: str, region: Optional[str], country: str):
        self.city: str = city
        self.state: Optional[str] = region
        self.country: str = country

    def __repr__(self) -> str:
        return f"Location({self.city}, {self.state}, {self.country})"


class Coordinates:
    """Coordinates represented by latitude and longitude (float values)"""

    def __init__(self, lat: float, long: float):
        self.lat: float = lat
        self.long: float = long

    def __repr__(self) -> str:
        return f"Coords({self.lat}, {self.long})"


class Paper:
    """Class representing paper"""

    def __init__(self, locations: list[Location]):
        self.locations: list[Location] = locations
