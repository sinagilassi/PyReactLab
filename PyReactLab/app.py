# import libs
from typing import List, Dict
# local
from pythermodb_settings.models import ComponentKey
from pyThermoLinkDB.models import ModelSource
from pyThermoLinkDB.thermo import Source
# locals
from .docs import (
    GasReactionSystem,
    LiquidReactionSystem
)
from .core import ReactionSystem, ReferenceManager


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
        model_source: ModelSource,
        component_key: ComponentKey = "Formula-State",
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
    model_source : ModelSource
        Inputs for the reaction system which contains the data source and equation source for the components in the reaction system.
    component_key : ComponentKey, optional
        The key to identify the components in the reaction system, by default "Formula-State".
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

        # NOTE: create source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # NOTE: create reaction system object
        return ReactionSystem(
            system_name=system_name,
            reactions=reactions,
            model_source=model_source,
            source=Source_
        )
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e


def create_gas_rxn(
    system_name: str,
    reactions: List[Dict[str, str]],
    model_source: ModelSource,
    component_key: ComponentKey = "Formula-State",
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
    model_source : ModelSource
        Inputs for the reaction system which
    component_key : ComponentKey, optional
        The key to identify the components in the reaction system, by default "Formula-State".
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

        # NOTE: create source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # NOTE: create reaction system object
        return GasReactionSystem(
            system_name=system_name,
            reactions=reactions,
            model_source=model_source,
            source=Source_
        )
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e


def create_liquid_rxn(
        system_name: str,
        reactions: List[Dict[str, str]],
        model_source: ModelSource,
        component_key: ComponentKey = "Formula-State",
        **kwargs
) -> ReactionSystem:
    """
    Create a liquid reaction system in that the compounds are in liquid phase and defined as `compound(l)` in the reaction expression.

    Parameters
    ----------
    system_name : str
        Name of the reaction system.
    reactions : list
        List of reactions in the system must be in the form of a list of dictionaries as the following keys:
        - 'reaction': str, the reaction equation.
        - 'name': str, the name of the reaction.
    model_source : ModelSource
        Inputs for the reaction system which
    component_key : ComponentKey, optional
        The key to identify the components in the reaction system, by default "Formula-State".
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

        # NOTE: create source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # NOTE: create reaction system object
        return LiquidReactionSystem(
            system_name=system_name,
            reactions=reactions,
            model_source=model_source,
            source=Source_
        )
    except Exception as e:
        raise Exception(f"Error creating reaction system: {e}") from e
