# import libs
from typing import List, Dict, Any, Literal, Optional
# local
from .docs import (
    ReactionSystem,
    Reaction,
    ReferenceManager,
    GasReactionSystem,
    LiquidReactionSystem
)
from .utils import model_source_checker
from .configs import DATASOURCE, EQUATIONSOURCE


def summary() -> str:
    """
    Provide a summary of the app calculation methods and required inputs.
    """
    try:
        return ReferenceManager().load_reference()
    except Exception as e:
        raise Exception(f"Error loading app summary: {e}") from e


def create_rxn(
        system_name: str,
        reactions: List[Dict[str, str]],
        model_source: Dict[str, Any],
        **kwargs
) -> ReactionSystem:
    """
    Create a reaction system.

    Parameters
    ----------
    system_name : str
        Name of the reaction system.
    reactions : list
        List of reactions in the system must be in the form of a list of dictionaries as the following keys:
        - 'reaction': str, the reaction equation.
        - 'name': str, the name of the reaction.
    model_source : dict
        Inputs for the reaction system which
    kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    ReactionSystem
        Reaction system object.
    """
    try:
        # NOTE: checking inputs
        # check if system name is valid
        if not isinstance(system_name, str):
            raise ValueError("System name must be a string.")

        # check if system inputs are valid
        if not isinstance(model_source, dict):
            raise ValueError("System inputs must be a dictionary.")

        # system source must contain datasource and equationsource
        # NOTE: check if model_source is valid
        if model_source:
            if not model_source_checker(model_source):
                raise ValueError(
                    f"Invalid model source, must contain {DATASOURCE} and {EQUATIONSOURCE}.")

        # NOTE: create reaction system object
        return ReactionSystem(
            system_name,
            reactions,
            model_source)
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e


def create_gas_rxn(
    system_name: str,
    reactions: List[Dict[str, str]],
    model_source: Dict[str, Any],
    **kwargs
) -> ReactionSystem:
    """
    Create a gaseous reaction system in that the compounds are in the gas phase and defined as `compound(g)`.

    Parameters
    ----------
    system_name : str
        Name of the reaction system.
    reactions : list
        List of reactions in the system must be in the form of a list of dictionaries as the following keys:
        - 'reaction': str, the reaction equation.
        - 'name': str, the name of the reaction.
    model_source : dict
        Inputs for the reaction system which
    kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    GasReactionSystem
        Gas Reaction system object.

    Notes
    -----
    The reaction system is created with the following attributes:

    - R1: `CO2(g)` + `3H2(g)` `<=>` `CH3OH(g)` + `H2O(g)`
    """
    try:
        # NOTE: checking inputs
        # check if system name is valid
        if not isinstance(system_name, str):
            raise ValueError("System name must be a string.")

        # check if system inputs are valid
        if not isinstance(model_source, dict):
            raise ValueError("System inputs must be a dictionary.")

        # system source must contain datasource and equationsource
        # NOTE: check if model_source is valid
        if model_source:
            if not model_source_checker(model_source):
                raise ValueError(
                    f"Invalid model source, must contain {DATASOURCE} and {EQUATIONSOURCE}.")

        # NOTE: create reaction system object
        return GasReactionSystem(
            system_name,
            reactions,
            model_source)
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e


def create_liquid_rxn(
        system_name: str,
        reactions: List[Dict[str, str]],
        model_source: Dict[str, Any],
        **kwargs
) -> ReactionSystem:
    """
    Create a liquid reaction system in that the compounds are in liquid phase and defined as `compound(l)`.

    Parameters
    ----------
    system_name : str
        Name of the reaction system.
    reactions : list
        List of reactions in the system must be in the form of a list of dictionaries as the following keys:
        - 'reaction': str, the reaction equation.
        - 'name': str, the name of the reaction.
    model_source : dict
        Inputs for the reaction system which
    kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    LiquidReactionSystem
        Liquid Reaction system object.

    Notes
    -----
    The reaction system is created with the following attributes:

    - R1: `CH3COOH(l)` + `C2H5OH(l)` `<=>` `CH3COOC2H5(l)` + `H2O(l)`
    """
    try:
        # NOTE: checking inputs
        # check if system name is valid
        if not isinstance(system_name, str):
            raise ValueError("System name must be a string.")

        # check if system inputs are valid
        if not isinstance(model_source, dict):
            raise ValueError("System inputs must be a dictionary.")

        # system source must contain datasource and equationsource
        # NOTE: check if model_source is valid
        if model_source:
            if not model_source_checker(model_source):
                raise ValueError(
                    f"Invalid model source, must contain {DATASOURCE} and {EQUATIONSOURCE}.")

        # NOTE: create reaction system object
        return LiquidReactionSystem(
            system_name,
            reactions,
            model_source)
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e
