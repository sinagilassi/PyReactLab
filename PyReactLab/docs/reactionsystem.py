# import libs
import numpy as np
from typing import Dict, Any, List, Literal, Optional
import pycuc
import time  # Import the time module
import pyThermoModels as ptm
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
    __reaction_list: Dict[str, Reaction] = {}

    # reference plugin
    _references = {}

    # overall reaction phase
    overall_reaction_phase = None

    def __init__(self,
                 system_name: str,
                 reactions: List[Dict[str, Any]],
                 model_source: Dict[str, Any],
                 **kwargs
                 ):
        self.__system_name = system_name
        self.__reactions = reactions

        # NOTE: model source
        self.__model_source = model_source

        # NOTE: kwargs
        # phase rule
        self.phase_rule = kwargs.get("phase_rule", None)

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
        # check
        if self.__system_name is None:
            return "No system name found."
        return self.__system_name

    @property
    def reactions(self) -> List[Dict[str, str]]:
        """Get the reactions of the reaction system."""
        # check
        if self.__reactions is None:
            raise ValueError("No reactions found.")
        return self.__reactions

    def select_reaction(self, reaction_name: str) -> Reaction:
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

            # reaction object
            reaction = self.__reaction_list.get(reaction_name, None)
            # check if reaction is valid
            if reaction is None or reaction == 'None':
                raise ValueError(
                    f"Invalid reaction object for {reaction_name}")
            return reaction
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
                r_ = Reaction(
                    self.datasource,
                    self.equationsource,
                    item,
                    phase_rule=self.phase_rule,)

                # NOTE: analyze reaction
                _res = r_.reaction_analysis_result

                # name
                name: str = item['name']
                # update
                reaction_res[name] = _res
                # set
                self.__reaction_list[name] = r_

            # NOTE: reaction numbers
            # self.reaction_numbers = len(reaction_res)

            # SECTION: analyze overall reaction
            # ! to set consumed, produced, and intermediate species
            # NOTE: version 1
            # res_0 = ChemReactUtils_.analyze_overall_reactions(
            #     self.reactions)
            # NOTE: version 2
            res_0 = ChemReactUtils_.analyze_overall_reactions_v2(
                reaction_res)

            # SECTION: set component
            # NOTE: version 1
            # res_1 = ChemReactUtils_.define_component_id(
            #     reaction_res)
            # NOTE: version 2
            res_1 = ChemReactUtils_.define_component_id_v2(
                reaction_res)

            # extract
            # ? component_list: list of components
            # ? component_dict: dict of component: id
            # ? comp_list: list of dict of component stoichiometry
            # ? comp_coeff: list of list of component stoichiometry
            component_list, component_dict, comp_list, comp_coeff = res_1

            # set stoichiometry transpose
            comp_coeff_t = np.array(comp_coeff).T

            # SECTION: overall reaction phase
            # NOTE: set reaction states
            overall_reaction_phases = []
            for item in reaction_res:
                # name
                name = item

                # NOTE: reaction state
                _res = self.__reaction_list[name].reaction_phase

                # save
                overall_reaction_phases.append(_res)

            # set
            overall_reaction_phases = list(set(overall_reaction_phases))

            # check if all reactions have the same phase
            if len(overall_reaction_phases) == 0:
                raise ValueError(
                    "No overall reaction phases found in the reactions.")

            # set
            if len(overall_reaction_phases) == 1:
                # set
                self.overall_reaction_phase = overall_reaction_phases[0]
            else:
                # set
                self.overall_reaction_phase = '-'.join(overall_reaction_phases)

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
            self.component_list = component_list
            self.component_dict = component_dict
            self.coeff_list_dict = comp_list
            self.coeff_list_list = comp_coeff
            self.coeff_T_list_list = comp_coeff_t
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
                                      **kwargs):
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
        **kwargs : dict
            Additional arguments for the calculation.
                - message: Optional message to display during the calculation.

        Returns
        -------
        res : dict
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
            res = reaction.cal_equilibrium_constant(
                temperature, method=method)

            # NOTE: message
            message = kwargs.get("message", None)
            if message is not None:
                res['message'] = message

            # res
            return res
        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.equilibrium_constant_at_temperature(): {str(e)}") from e

    def equilibrium(self,
                    inputs: Dict[str, Any],
                    conversion: Optional[List[str]] = None,
                    gas_mixture: Literal[
                        "ideal", "non-ideal"
                    ] = "ideal",
                    solution: Literal[
                        "ideal", "non-ideal"
                    ] = "ideal",
                    method: Literal[
                        'minimize', 'least_squares'
                    ] = 'minimize',
                    **kwargs):
        """
        Calculate the equilibrium state of the reaction system.

        Parameters
        ----------
        inputs : dict
            Inputs for the equilibrium calculation.
        conversion : list, optional
            List of components to calculate conversion for, by default None.
        gas_mixture : str, optional
            Type of gas mixture, by default "ideal".
        solution : str, optional
            Type of liquid mixture, by default "ideal".
        method : str, optional
            Method for the calculation, by default "minimize".
            Options are "minimize" or "least_squares".
        **kwargs : dict
            Additional arguments for the calculation.
                - eos_model: Equation of state model to use for the calculation. Options are "SRK" or "PR".
                - activity_model: Activity model to use for the calculation. Options are "NRTL" or "UNIFAC".
                - message: Optional message to display during the calculation.

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

        For `non-ideal solution`, the rules of the activity model should be followed.

        Activity inputs for NRTL model should contain:
        - tau: dict, tau values for the components.
        - alpha: dict, alpha values for the components.

        Activity inputs for UNIQUAC model should contain:
        - q: dict, q values for the components.
        - r: dict, r values for the components.
        - tau: dict, tau values for the components.

        For both models, dU or a, b, c, and d are required to calculate tau in case tau is not provided. All values should be introduced in the `inputs` as:

        ```python
        inputs = {
            'mole': mole,
            'temperature': [100, "C"],
            'pressure': [1.0, "bar"],
            'tau': tau,
            'alpha': alpha, ...
            }
        ```
        """
        try:
            # ! Start timing
            start_time = time.time()

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

            # NOTE: build mole and mole fraction matrix regarding the component dict
            # check
            if initial_mole is not None and initial_mole_fraction is not None:
                # set values
                initial_mole_std, initial_mole_fraction_std, _, _ = ReactionAnalyzer.set_stream(
                    component_dict=self.component_dict,
                    mole=initial_mole,
                    mole_fraction=initial_mole_fraction,
                )
            else:
                raise ValueError(
                    "Initial mole and mole fraction must be provided.")

            # SECTION: kwargs
            # eos model (name)
            eos_model = kwargs.get("eos_model", "SRK")
            # check if eos model is valid
            if eos_model not in ["SRK", "PR", 'RK', 'vdW']:
                raise ValueError(
                    "Invalid eos model. Options are 'SRK','RK','PR', or 'vdW'.")

            # activity model (name)
            activity_model = kwargs.get("activity_model", "NRTL")
            # check if activity model is valid
            if activity_model not in ["NRTL", "UNIQUAC"]:
                raise ValueError(
                    "Invalid activity model. Options are 'NRTL' or 'UNIQUAC'.")

            # SECTION: init
            ReactionOptimizer_ = ReactionOptimizer(
                self.datasource,
                self.equationsource,
                self.component_dict,
                self.coeff_list_dict,
                self.coeff_T_list_list,
                self.reaction_analysis,
                self.overall_reaction_analysis,
            )

            # SECTION: setting up the reaction optimizer
            # NOTE: set up eos model
            # ? eos model
            ReactionOptimizer_.eos_model = eos_model
            # init eos class
            ReactionOptimizer_.eos = ptm.eos()

            # gas mixture
            ReactionOptimizer_.gas_mixture = gas_mixture

            # NOTE: set up activity model
            # ? activity inputs
            # NRTL: tau, alpha,
            # dg or a, b, c, and d are required to calculated tau
            # UNIQUAC: q, r, tau
            # dU or a, b, c, and d are required to calculated tau
            activity_inputs = {**inputs}
            # remove temperature and pressure
            activity_inputs.pop("temperature", None)
            activity_inputs.pop("pressure", None)
            # remove mole and mole fraction
            activity_inputs.pop("mole", None)
            activity_inputs.pop("mole_fraction", None)

            # set
            if activity_inputs is not None:
                #  check if activity inputs is valid
                if not isinstance(activity_inputs, dict):
                    raise ValueError(
                        "Activity inputs must be a dictionary.")

                # set
                ReactionOptimizer_.activity_inputs = activity_inputs

            # ? activity model
            ReactionOptimizer_.activity_model = activity_model
            # init activity class
            ReactionOptimizer_.activity = ptm.activities(
                components=self.component_list,
                model_name=activity_model,
                model_source=self.model_source,
            )

            # solution
            ReactionOptimizer_.solution = solution

            # SECTION: equilibrium constant calculation
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
            opt_res = ReactionOptimizer_.opt_run(
                initial_mole=initial_mole_std,
                initial_mole_fraction=initial_mole_fraction_std,
                temperature=temperature_K,
                pressure=pressure_bar,
                equilibrium_constants=equilibrium_constant,
                reaction_numbers=self.reaction_numbers,
                method=method,
            )

            # NOTE: check optimization result
            res = {}

            # check if optimization is successful
            if opt_res.success:
                # extent of reaction
                EoR = opt_res.x

                # calculate equilibrium state
                Eq_Xfs, Eq_Xfs_vector, Eq_Nf, Eq_Nfs_vector, Eq_Nfs = \
                    ReactionOptimizer_.equilibrium_results(
                        initial_mole=initial_mole_std,
                        EoR=EoR,
                    )

                # NOTE: calculate conversion
                # res
                conversion_res = None
                # conversion
                if conversion is not None:
                    # check list
                    if not isinstance(conversion, list):
                        raise ValueError(
                            "Conversion must be a list of strings.")

                    # check component in the components
                    for key in conversion:
                        if key not in self.component_dict:
                            raise ValueError(
                                f"Invalid component: {key} in the components.")

                    # calculate conversion
                    conversion_res = ReactionAnalyzer.cal_conversion(
                        initial_mole=initial_mole_std,
                        final_mole=Eq_Nfs,
                        components=conversion,
                    )

                # SECTION: set results
                # initial feed
                res['feed'] = {
                    "mole": {
                        'value': initial_mole_std,
                        'unit': "mol"
                    },
                    "mole_total": {
                        'value': sum(initial_mole_std.values()),
                        'unit': "mol"
                    },
                    "mole_fraction": {
                        'value': initial_mole_fraction_std,
                        'unit': "dimensionless"
                    },
                    "mole_fraction_sum": {
                        'value': sum(initial_mole_fraction_std.values()),
                        'unit': "dimensionless"
                    },
                }

                # equilibrium state
                res['equilibrium'] = {
                    "mole": {
                        'value': Eq_Nfs,
                        'unit': "mol"
                    },
                    "mole_total": {
                        'value': sum(Eq_Nfs.values()),
                        'unit': "mol"
                    },
                    "mole_fraction": {
                        'value': Eq_Xfs,
                        'unit': "dimensionless"
                    },
                    "mole_fraction_sum": {
                        'value': sum(Eq_Xfs.values()),
                        'unit': "dimensionless"
                    },
                }

                # extent of reaction
                res['extent_of_reaction'] = {
                    "value": EoR,
                    "unit": "dimensionless"
                }

                # equilibrium condition
                # temperature
                res['temperature'] = {
                    "value": temperature_K,
                    "unit": "K"
                }

                # pressure
                res['pressure'] = {
                    "value": pressure_bar,
                    "unit": "bar"
                }

                # optimization fun value
                res['optimization_fun'] = {
                    "value": opt_res.fun,
                    "unit": "dimensionless"
                }

                # conversion
                if conversion_res is not None:
                    res['conversion'] = conversion_res

                # message
                message = kwargs.get("message", None)
                if message is not None:
                    res['message'] = message

                # gas mixture
                res['gas_mixture'] = gas_mixture
                # solution
                res['solution'] = solution

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # res
            return res

        except Exception as e:
            raise Exception(
                f"Failing in the equilibrium calculations {str(e)}") from e
