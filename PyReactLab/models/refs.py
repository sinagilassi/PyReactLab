# create reference
from typing import Dict, TypedDict
from pythermodb_settings.models import (
    Temperature,
    Pressure,
)


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
