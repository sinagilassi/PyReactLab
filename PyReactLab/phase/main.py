# import libs
from typing import List, Dict, Any
import logging
# local
from ..models import Temperature, Pressure, Component
from .phasecontroller import PhaseController


# NOTE: get logger
logger = logging.getLogger(__name__)


def check_phase(
    component: Component,
    temperature: Temperature,
    pressure: Pressure,
    datasource: Dict[str, Any],
    equationsource: Dict[str, Any],
):
    """
    Check the phase of a single component.

    Parameters
    ----------
    component : Component
        The component object containing name, formula, state, and mole_fraction.
    temperature : Temperature
        The temperature object containing value and unit as [300, 'K'].
    pressure : Pressure
        The pressure object containing value and unit as [101325, 'Pa'].
    datasource : Dict[str, Any]
        The datasource dictionary containing thermodynamic data.
    equationsource : Dict[str, Any]
        The equationsource dictionary containing equations for calculations.

    Returns
    -------
    str
        The phase of the component: 'solid', 'liquid', or 'gas'.
    """
    try:
        # create phase controller
        phase_controller = PhaseController(
            datasource=datasource,
            equationsource=equationsource
        )
        # check phase
        phase = phase_controller.check_component_phase(
            temperature=temperature,
            pressure=pressure,
            component=component
        )
        return phase
    except Exception as e:
        logger.error(f"Error in check_phase: {e}")
        raise
