# import libs
from typing import Dict, Any, List, Literal, Optional
import pycuc
# local
from .reaction import Reaction
from .thermolinkdb import ThermoLinkDB
from .refmanager import ReferenceManager
from .reactionanalyzer import ReactionAnalyzer
from .optim import ReactionOptimizer
from ..utils import ChemReactUtils
from .reaction import Reaction


class ReactionSystem(ThermoLinkDB, ReferenceManager):
    """Class to represent a system of chemical reactions."""

    # NOTE: class variables
    __system_name = None
    __reactions = None
    # primary analysis result
    __reaction_analysis = None
    __reaction_list = {}

    # reference plugin
    _references = {}

    def __init__(self,
                 system_name: str,
                 reactions: List[Dict[str, Any]],
                 model_source: Dict[str, Any]
                 ):
        self.__system_name = system_name
        self.__reactions = reactions

        # NOTE: model source
        self.__model_source = model_source

        # NOTE: init class
        ReferenceManager.__init__(self)
        ThermoLinkDB.__init__(self, model_source)

        # SECTION: load reference
        # reference plugin (default app params)
        self._references = self.load_reference()

        # SECTION: energy analysis result list
        self.__reaction_analyzer()

    @property
    def system_name(self) -> str:
        """Get the name of the reaction system."""
        return self.__system_name

    @property
    def reactions(self) -> List[Dict[str, Any]]:
        """Get the reactions of the reaction system."""
        return self.__reactions

    def select_reaction(self, reaction_name: str) -> Optional[Reaction]:
        """
        Select a reaction from the reaction system.

        Parameters
        ----------
        reaction_name : str
            Name of the reaction.

        Returns
        -------
        Reaction or None
            Reaction object if found, otherwise None.
        """
        try:
            # NOTE: check if reaction name is valid
            if reaction_name not in self.__reaction_list:
                raise ValueError(f"Invalid reaction name: {reaction_name}")

            # return reaction object
            return self.__reaction_list[reaction_name]
        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.select_reaction(): {str(e)}") from e

    def __reaction_analyzer(self) -> None:
        """
        Execute the primary analysis for the reaction system.
        """
        try:
            # NOTE: initialize
            ChemReactUtils_ = ChemReactUtils()
            # ReactionAnalyzer_ = ReactionAnalyzer()

            # SECTION: reaction system analysis
            # analyze reaction
            reaction_res = {}

            # looping through each reaction
            for item in self.reactions:
                # NOTE: create reaction
                r_ = Reaction(self.datasource, self.equationsource, item)

                # NOTE: analyze reaction
                _res = r_.reaction_analysis_result

                # name
                name = item['name']
                # update
                reaction_res[name] = _res
                # set
                self.__reaction_list[name] = r_

            # NOTE: reaction numbers
            # self.reaction_numbers = len(reaction_res)

            # SECTION: analyze overall reaction
            # ! to set consumed, produced, and intermediate species
            res_0 = ChemReactUtils_.analyze_overall_reactions(
                self.reactions)

            # SECTION: set component
            res_1 = ChemReactUtils_.define_component_id(
                reaction_res)
            # extract
            component_list, component_dict, comp_list, comp_coeff = res_1

            # SECTION: energy analysis
            # energy analysis result list
            energy_analysis = {}

            # loop through each reaction
            for item in reaction_res:
                # name
                name = item

                # NOTE: energy analysis
                reaction_ = self.__reaction_list[name]
                _res = reaction_.energy_analysis

                # save
                energy_analysis[item] = _res

            # NOTE: set primary analysis result
            self.reaction_numbers = len(reaction_res)
            self.reaction_analysis = reaction_res
            self.overall_reaction_analysis = res_0
            self.reaction_states = None
            self.component_list = component_list  # ? list of components
            self.component_dict = component_dict  # ? dict of component: id
            self.coeff_list_dict = comp_list  # ? list of dict of component: coeff
            self.coeff_list_list = comp_coeff  # ? list of list of component coeff
            self.energy_analysis = energy_analysis

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.go(): {str(e)}") from e

    def reaction_equilibrium_constant(self,
                                      reaction_name: str,
                                      temperature: list[float | str],
                                      method: Literal[
                                          "van't Hoff", "shortcut van't Hoff"
                                      ] = "van't Hoff",
                                      ):
        """
        Calculate the equilibrium constant at a given temperature using the van't Hoff equation.

        Parameters
        ----------
        reaction_name : str
            Name of the reaction.
        temperature : list
            Temperature in the form of [value, unit], the unit is automatically converted to K.
        method : str, optional
            Method to calculate the equilibrium constant, by default "van't Hoff".
            Options are "van't Hoff" or "shortcut van't Hoff".

        Returns
        -------
        float
            Equilibrium constant at the given temperature.
        """
        try:
            # NOTE: check if temperature is valid
            if not isinstance(temperature, list):
                raise ValueError("Temperature must be a number.")

            # check if temperature is valid
            if len(temperature) != 2:
                raise ValueError(
                    "Temperature must be a list of length 2 following [value, unit].")

            if not isinstance(temperature[0], (int, float)):
                raise ValueError("Temperature must be a number.")

            if not isinstance(temperature[1], str):
                raise ValueError("Temperature unit must be a string.")

            # NOTE: check if reaction name is valid
            if reaction_name not in self.__reaction_list:
                raise ValueError(f"Invalid reaction name: {reaction_name}")

            # SECTION: pre-calculation
            # get reaction object
            reaction = self.__reaction_list[reaction_name]

            # check if reaction is valid
            if not isinstance(reaction, Reaction):
                raise ValueError(
                    f"Invalid reaction object for {reaction_name}")

            # NOTE: calculate equilibrium constant at the given temperature
            return reaction.cal_equilibrium_constant(
                temperature, method=method)

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.equilibrium_constant_at_temperature(): {str(e)}") from e

    def equilibrium(self,
                    inputs: Dict[str, Any],
                    gas_mixture: Literal[
                        "ideal", "non-ideal"
                    ] = "ideal",
                    liquid_mixture: Literal[
                        "ideal", "non-ideal"
                    ] = "ideal",
                    **kwargs):
        """
        Calculate the equilibrium state of the reaction system.

        Parameters
        ----------
        inputs : dict
            Inputs for the equilibrium calculation which must contain:
        **kwargs : dict
            Additional arguments for the calculation.
                - eos_model: Equation of state model to use for the calculation. Options are "SRK" or "PR".
                - activity_model: Activity model to use for the calculation. Options are "NRTL" or "UNIFAC".

        Returns
        -------
        dict
            Equilibrium state of the reaction system.

        Notes
        -----
        Inputs must contain `mole_fraction` or `mole`, `temperature`, and `pressure` as follows:

        - mole_fraction: dict, initial mole fraction of the components in the system.
        - mole: dict, initial mole of the components in the system.
        - temperature: list, temperature in the form of [value, unit].
        - pressure: list, pressure in the form of [value, unit].
        """
        try:
            # SECTION: check args
            # NOTE: check if inputs are valid
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")
            # check if inputs are valid
            if not any(key in inputs for key in ["mole_fraction", "mole"]):
                raise ValueError(
                    "Inputs must contain either mole_fraction or mole.")
            if "temperature" not in inputs:
                raise ValueError("Inputs must contain temperature.")
            if "pressure" not in inputs:
                raise ValueError("Inputs must contain pressure.")

            # NOTE: check if temperature is valid
            # set
            temperature = inputs["temperature"]
            # check if temperature is valid
            if not isinstance(temperature, list):
                raise ValueError("Temperature must be a number.")

            # check if temperature is valid
            if len(temperature) != 2:
                raise ValueError(
                    "Temperature must be a list of length 2 following [value, unit].")

            if not isinstance(temperature[0], (int, float)):
                raise ValueError("Temperature must be a number.")

            if not isinstance(temperature[1], str):
                raise ValueError("Temperature unit must be a string.")

            # ! convert to K
            # set unit
            unit_set = f"{temperature[1]} => K"
            temperature_K = pycuc.to(temperature[0], unit_set)

            # NOTE: check if pressure is valid
            # set
            pressure = inputs["pressure"]
            # check if pressure is valid
            if not isinstance(pressure, list):
                raise ValueError("Pressure must be a number.")

            # check if pressure is valid
            if len(pressure) != 2:
                raise ValueError(
                    "Pressure must be a list of length 2 following [value, unit].")

            if not isinstance(pressure[0], (int, float)):
                raise ValueError("Pressure must be a number.")

            if not isinstance(pressure[1], str):
                raise ValueError("Pressure unit must be a string.")

            # ! convert to bar
            # set unit
            unit_set = f"{pressure[1]} => bar"
            pressure_bar = pycuc.to(pressure[0], unit_set)

            # NOTE: check if initial mole fraction is valid
            # set
            initial_mole_fraction = inputs.get("mole_fraction", None)
            # check
            initial_mole = inputs.get("mole", None)

            if initial_mole_fraction is not None:
                # check if initial mole fraction is valid
                if not isinstance(initial_mole_fraction, dict):
                    raise ValueError(
                        "Initial mole fraction must be a dictionary.")

                # check if initial mole fraction is valid
                for key in initial_mole_fraction:
                    if not isinstance(key, str):
                        raise ValueError(
                            "Initial mole fraction key must be a string.")
                    if not isinstance(initial_mole_fraction[key], (int, float)):
                        raise ValueError(
                            "Initial mole fraction value must be a number.")

                # ? normalize
                # check if sum of mole fraction is 1
                initial_mole_fraction = ReactionAnalyzer.norm_mole_fraction(
                    initial_mole_fraction)

                # NOTE: convert to mole fraction
                initial_mole = ReactionAnalyzer.cal_mole(
                    initial_mole_fraction)

            elif initial_mole is not None:
                # check if initial mole is valid
                if not isinstance(initial_mole, dict):
                    raise ValueError("Initial mole must be a dictionary.")

                # check if initial mole is valid
                for key in initial_mole:
                    if not isinstance(key, str):
                        raise ValueError(
                            "Initial mole key must be a string.")
                    if not isinstance(initial_mole[key], (int, float)):
                        raise ValueError(
                            "Initial mole value must be a number.")

                # NOTE: convert to mole fraction
                initial_mole_fraction, _ = ReactionAnalyzer.cal_mole_fraction(
                    initial_mole)

            # SECTION: kwargs
            # eos model
            eos_model = kwargs.get("eos_model", "SRK")

            # activity model
            activity_model = kwargs.get("activity_model", "NRTL")

            # SECTION: init
            ReactionOptimizer_ = ReactionOptimizer(
                self.datasource,
                self.equationsource,
                self.component_dict,
                self.coeff_list_dict,
                self.reaction_analysis,
                self.overall_reaction_analysis,
            )

            # NOTE: setting up the reaction optimizer
            # eos model
            ReactionOptimizer_.eos_model = eos_model
            # init eos class
            ReactionOptimizer_.init_eos()
            # activity model
            ReactionOptimizer_.activity_model = activity_model
            # init activity class
            ReactionOptimizer_.init_activity()

            # NOTE: equilibrium constant calculation
            # res
            equilibrium_constant = {}

            # loop through each reaction
            for key, r in self.__reaction_list.items():
                # NOTE: check if reaction is valid
                if not isinstance(r, Reaction):
                    raise ValueError(
                        f"Invalid reaction object for {key}")

                # ! calculate equilibrium constant at the given temperature
                res_ = r.cal_equilibrium_constant(temperature)

                # save
                equilibrium_constant[key] = res_

            # NOTE: run equilibrium calculation
            res = ReactionOptimizer_.opt_run(
                initial_mole=initial_mole,
                temperature=temperature_K,
                pressure=pressure_bar,
                equilibrium_constants=equilibrium_constant,
                reaction_numbers=self.reaction_numbers,
            )

            # res
            return res

        except Exception as e:
            raise Exception(
                f"Failing in the equilibrium calculations {str(e)}") from e
