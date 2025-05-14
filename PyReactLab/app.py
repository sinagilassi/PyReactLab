# import libs
from typing import List, Dict, Any, Literal, Optional
# local
from .docs import ReactionSystem, Reaction, ReferenceManager
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
        reactions: List[Dict[str, Any]],
        model_source: Dict[str, Any],
        **kwargs) -> ReactionSystem:
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
