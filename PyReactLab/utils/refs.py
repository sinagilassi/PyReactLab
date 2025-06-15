# create reference
from typing import Dict, TypedDict


class Temperature(TypedDict):
    """
    Temperature reference for the reaction system.
    """

    value: float  # Temperature in Kelvin
    unit: str  # Unit of temperature, e.g., 'K', 'C'


class Pressure(TypedDict):
    """
    Pressure reference for the reaction system.
    """

    value: float  # Pressure in Pascal
    unit: str  # Unit of pressure, e.g., 'Pa', 'bar', 'atm'


class OperatingConditions(TypedDict):
    """
    Operating conditions for the reaction system.
    """

    temperature: Temperature  # Temperature in Kelvin
    pressure: Pressure  # Pressure in Pascal


class MoleFraction(TypedDict):
    """
    Mole fraction reference for the reaction system.
    """

    # Dictionary of species and their mole fractions
    components: Dict[str, float | int]
