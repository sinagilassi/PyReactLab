# import libs
import logging
from typing import Dict, Any, List
import pycuc
# local
from ..models import Temperature, Pressure, Component
from ..configs import (
    DATASOURCE,
    EQUATIONSOURCE,
    BOILING_TEMPERATURE,
    MELTING_TEMPERATURE,
    SUBLIMATION_TEMPERATURE
)
from ..thermo import ThermoSource

# NOTE: get logger
logger = logging.getLogger(__name__)


class PhaseController:
    """
    PhaseController class for managing phase-related operations in PyReactLab.

    This class provides methods to handle phase calculations, the datasource and equationsource must include the necessary thermodynamic data such as:
    - Boiling temperature (Tb)
    - Melting temperature (Tm)
    - Sublimation temperature (Ts)
    - Heat of vaporization (EnVap)
    - Heat of fusion (EnFus)
    - Heat of sublimation (EnSub)
    - Heat capacity at constant pressure for gas (Cp_GAS), liquid (Cp_LIQ), and solid (Cp_SOL)
    """

    def __init__(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any]
    ):
        """
        Initialize the PhaseController instance.
        """
        # LINK to ThermoSource
        self.ThermoSource_ = ThermoSource(model_source={
            DATASOURCE: datasource,
            EQUATIONSOURCE: equationsource
        })

    def check_component_phase(
        self,
        temperature: Temperature,
        pressure: Pressure,
        component: Component,
        component_mode: str = 'formula'
    ) -> str:
        """
        Check the component phase.

        Parameters
        ----------
        temperature : Temperature
            The temperature object containing value and unit as [298.15, 'K'].
        pressure : Pressure
            The pressure object containing value and unit as [101325, 'Pa'].
        component : Component
            The component object containing name, formula, state, and mole_fraction.
        component_mode : str, optional
            The mode to identify the component, either 'name' or 'formula'. Default is 'formula'.

        Returns
        -------
        str
            The phase of the component: 'g' for gas, 'l' for liquid, 's' for solid.
        """
        try:
            # SECTION: extract properties
            # ! component id
            component_id = component.formula if component_mode == 'formula' else component.name
            # add state
            component_id = f"{component_id}-{component.state}"

            # ! boiling temperature
            Tb_res = self.ThermoSource_.data_extractor(
                component_name=component_id,
                prop_name=BOILING_TEMPERATURE
            )

            # get value and unit
            Tb_value = Tb_res.get('value', None)
            Tb_unit = Tb_res.get('unit', None)

            # check value and unit
            if (
                Tb_value is not None and
                Tb_unit is not None
            ):
                # unit conversion
                Tb_value = pycuc.to(Tb_value, f"{Tb_unit} => K")

                # update
                Tb = Temperature(
                    value=Tb_value,
                    unit="K"
                )
            else:
                raise ValueError("Boiling temperature data is incomplete.")

            # ! melting temperature
            Tm_res = self.ThermoSource_.data_extractor(
                component_name=component_id,
                prop_name=MELTING_TEMPERATURE
            )

            # get value and unit
            Tm_value = Tm_res.get('value', None)
            Tm_unit = Tm_res.get('unit', None)

            # check value and unit
            if (
                Tm_value is not None and
                Tm_unit is not None
            ):
                # unit conversion
                Tm_value = pycuc.to(Tm_value, f"{Tm_unit} => K")

                # update
                Tm = Temperature(
                    value=Tm_value,
                    unit="K"
                )
            else:
                raise ValueError("Melting temperature data is incomplete.")

            # ! sublimation temperature
            Ts_res = self.ThermoSource_.data_extractor(
                component_name=component_id,
                prop_name=SUBLIMATION_TEMPERATURE
            )

            # get value and unit
            Ts_value = Ts_res.get('value', None)
            Ts_unit = Ts_res.get('unit', None)

            # check value and unit
            if (
                Ts_value is not None and
                Ts_unit is not None
            ):
                # unit conversion
                Ts_value = pycuc.to(Ts_value, f"{Ts_unit} => K")

                # update
                Ts = Temperature(
                    value=Ts_value,
                    unit="K"
                )
            else:
                raise ValueError("Sublimation temperature data is incomplete.")

            # SECTION: current temperature and pressure
            # unit conversion
            # kelvin [K]
            T = pycuc.to(temperature.value, f"{temperature.unit} => K")
            # pascal [Pa]
            P = pycuc.to(pressure.value, f"{pressure.unit} => Pa")

            # SECTION: check phase
            if T >= Tb.value:
                return 'g'  # gas
            elif Tm.value <= T < Tb.value:
                return 'l'  # liquid
            elif Ts.value <= T < Tm.value:
                return 's'  # solid
            else:
                raise ValueError("Temperature is below sublimation point.")

        except Exception as e:
            logger.error(f"Error checking component phase: {e}")
            raise

    def check_components_phase(
        self,
        temperature: Temperature,
        pressure: Pressure,
        components: List[Component],
        component_mode: str = 'formula'
    ) -> Dict[str, str]:
        """
        Check the phases of multiple components.

        Parameters
        ----------
        temperature : Temperature
            The temperature object containing value and unit as [298.15, 'K'].
        pressure : Pressure
            The pressure object containing value and unit as [101325, 'Pa'].
        components : List[Component]
            A list of component objects each containing name, formula, state, and mole_fraction.
        component_mode : str, optional
            The mode to identify the components, either 'name' or 'formula'. Default is 'formula'.

        Returns
        -------
        Dict[str, str]
            A dictionary mapping component names to their phases: 'g' for gas, 'l' for liquid, 's' for solid.
        """
        try:
            phases = {}
            for comp in components:
                phase = self.check_component_phase(
                    temperature=temperature,
                    pressure=pressure,
                    component=comp,
                    component_mode=component_mode
                )
                phases[comp.name] = phase
            return phases
        except Exception as e:
            logger.error(f"Error checking components phase: {e}")
            raise
