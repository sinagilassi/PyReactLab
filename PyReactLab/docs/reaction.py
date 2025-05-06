# import libs
from typing import List, Dict, Any, Literal, Optional
# local


class Reaction():
    """Class to represent a chemical reaction."""

    # variables
    __reactants = None
    __products = None

    def __init__(self, reactants, products, conditions=None):
        self.__reactants = reactants
        self.__products = products

    def __str__(self):
        """String representation of the reaction."""
        reactants_str = " + ".join(self.__reactants)
        products_str = " + ".join(self.__products)
        return f"{reactants_str} -> {products_str}]"

    def __repr__(self):
        """String representation of the reaction."""
        return self.__str__()

    @property
    def reactants(self):
        """Get the reactants of the reaction."""
        return self.__reactants

    @property
    def products(self):
        """Get the products of the reaction."""
        return self.__products

    @property
    def conditions(self):
        """Get the conditions of the reaction."""
        return self.__conditions
