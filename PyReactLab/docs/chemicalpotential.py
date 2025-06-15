# import libs
import numpy as np
from typing import Literal, Dict, Any, List
from math import log
# local
# local
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K, EOS_MODELS, ACTIVITY_MODELS
)
from .reaction import Reaction
from .reactionanalyzer import ReactionAnalyzer
from ..utils import (
    Temperature,
    Pressure,
    OperatingConditions,
)


class ChemicalPotential:
    """
    Class to calculate chemical potential of a reaction system.

    Notes
    -----
    The fugacity of an ideal-mixture gas component is defined as:

        f_i_ID = y_i * f_i_ID_PURE = y_i * P

    The Fugacity of an ideal-solution liquid component is defined as:

        f_i_ID(T,P,x) = x_i * f_i(T,P)

    f_i(T,P) is calculated using EOS (the smallest Z) or Poynting equation.

    The chemical potential of ideal-solution liquid component is defined as:

        μ_i_ID(T,P,x) = GiEnFo(T,P) + RT * ln(f_i_ID(T,P,x)/f_i(T,P))
        μ_i_ID(T,P,x) = GiEnFo(T,P) + RT * ln(x_i)

    The chemical potential of ideal-gas component is defined as:

        μ_i_ID(T,P,y) = GiEnFo(T,P) + RT * ln(f_i_ID(T,P,y)/f_i(T,P))
        μ_i_ID(T,P,y) = GiEnFo(T,P) + RT * ln(y_i)

    The chemical potential of non-ideal-gas component is defined as:

        μ_i_NID(T,P,y) = μ_i_ID(T,P,y) + RT * ln(φ_i(T,P,y))

        φ_i(T,P,y) = f_i_NID(T,P,y) / f_i_ID(T,P,y)

        μ_i_NID(T,P,y) = GiEnFo(T,P) + RT * ln(y_i) + RT * ln(φ_i(T,P,y))

    The chemical potential of non-ideal-solution component is defined as:

        μ_i_NID(T,P,x) = μ_i_ID(T,P,x) + RT * ln(AcCo_i(T,P,x))

        AcCo_i(T,P,x) = f_i_NID(T,P,x) / f_i_ID(T,P,x)

        μ_i_NID(T,P,x) = GiEnFo(T,P) + RT * ln(x_i) + RT * ln(AcCo_i(T,P,x))
    """
    # SECTION: attributes
    # system name
    _system_inputs = None
    # universal gas constant [J/mol.K]
    R = R_CONST_J__molK
    # temperature [K]
    T_Ref_K = TEMPERATURE_REF_K
    # pressure [Pa]
    P_Ref_Pa = PRESSURE_REF_Pa
    # pressure [bar]
    P_Ref_bar = 1

    # NOTE: models
    # eos
    _eos = None
    # eos model
    _eos_model = None
    # activity model
    _activity_model = None
    # activity
    _activity = None
    # activity inputs
    _activity_inputs = None

    # NOTE: systems
    # gas phase
    _gas_mixture = 'ideal'
    # solution
    _solution = 'ideal'

    def __init__(self,
                 datasource: Dict[str, Any],
                 equationsource: Dict[str, Any],
                 reaction_list: Dict[str, Reaction],
                 component_dict: Dict[str, float | int],
                 component_state_list: List[tuple],
                 comp_list: List[Dict[str, float | int]],
                 reaction_analysis: Dict,
                 phase_stream: Dict[str, Any],
                 phase_contents: Dict[str, Any],
                 overall_reaction_analysis: Dict,
                 overall_reaction_phase: str,
                 **kwargs):
        '''
        Initialize ReactionOptimizer class

        Parameters
        ----------
        datasource : dict
            datasource dictionary
        equationsource : dict
            equationsource dictionary
        component_dict : dict
            component dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        comp_list : list
            component list [{key: value}, {key: value}, ...], such as [{CO2: -1.0, H2: -3.0, CH3OH: 1.0, H2O: 1.0}, ...]
        reaction_analysis : dict
            reaction analysis result
        overall_reaction_analysis : dict
            overall reaction analysis result
        overall_reaction_phase : str
            overall reaction phase, such as 'gas', 'liquid', 'solid', or 'aqueous'
        kwargs : dict
            additional parameters
        '''
        # datasource
        self.datasource = datasource
        # equationsource
        self.equationsource = equationsource
        # set component dictionary
        self.component_dict = component_dict
        # set component state list
        self.component_state_list = component_state_list
        # reaction list
        self.reaction_list = reaction_list
        # set component list
        self.comp_list = comp_list
        # set reaction analysis
        self.reaction_analysis = reaction_analysis
        # phase stream
        self.phase_stream = phase_stream
        # set phase contents
        self.phase_contents = phase_contents
        # set overall reaction analysis
        self.overall_reaction_analysis = overall_reaction_analysis
        # set overall reaction phase
        self.overall_reaction_phase = overall_reaction_phase

        # NOTE: set
        # component list
        self.component_list = list(component_dict.keys())
        # component index
        self.component_index = list(component_dict.values())

        # NOTE: kwargs
        self.threshold = kwargs.get('threshold', 1e-8)

        # SECTION: init class
        self.ReactionAnalyzer_ = ReactionAnalyzer()

    @property
    def eos(self):
        '''Get eos class'''
        # check
        return self._eos

    @eos.setter
    def eos(self, value: Any):
        '''Set eos class'''
        # set
        self._eos = value

    @property
    def eos_model(self):
        '''Get eos model'''
        if self._eos_model is None:
            # set default eos model
            self._eos_model = 'SRK'
        return self._eos_model

    @eos_model.setter
    def eos_model(self, value: str):
        '''Set eos model'''
        # set
        self._eos_model = value

    @property
    def activity_model(self):
        '''Get activity model name'''
        if self._activity_model is None:
            # set default activity model
            self._activity_model = 'NRTL'
        return self._activity_model

    @activity_model.setter
    def activity_model(self, value: str):
        '''Set activity model name'''
        # set
        self._activity_model = value

    @property
    def activity(self):
        '''activity model (NRTL, UNIQUAC)'''
        return self._activity

    @activity.setter
    def activity(self, value: Any):
        '''Set activity model'''
        # set
        self._activity = value

    @property
    def activity_inputs(self):
        '''Get activity inputs'''
        return self._activity_inputs

    @activity_inputs.setter
    def activity_inputs(self, value: Dict[str, Any]):
        '''Set activity inputs'''
        # set
        self._activity_inputs = value

    @property
    def gas_mixture(self):
        '''Get gas mixture'''
        return self._gas_mixture

    @gas_mixture.setter
    def gas_mixture(self, value: str):
        '''Set gas mixture'''
        # set
        self._gas_mixture = value

    @property
    def solution(self):
        '''Get solution'''
        return self._solution

    @solution.setter
    def solution(self, value: str):
        '''Set solution'''
        # set
        self._solution = value

    def _cal_fugacity_coefficient_gaseous_mixture(
        self,
        model_name: str,
        model_input: Dict[str, Any]
    ):
        """
        Calculate the fugacity coefficient of gaseous mixture using the specified EOS model.

        Parameters
        ----------
        model_name : str
            The EOS model to use for calculation. Options: "SRK", "PR", "RK".
        model_input : dict
            The input data for the EOS model.
            - feed-specification: Dictionary of component mole fractions.
            - pressure: pressure value in any unit.
            - temperature: temperature value in any unit.

        Returns
        -------
        dict
            A dictionary containing the fugacity coefficients for each component in the gaseous mixture.
        """
        try:
            # SECTION: model source
            model_source = {
                "datasource": self.datasource,
                "equationsource": self.equationsource
            }

            # NOTE: calculate fugacity
            eos = self.eos

            # check
            if eos is None or eos == 'None':
                raise ValueError(
                    f"Invalid EOS model. Must be {EOS_MODELS}.")

            res = eos.cal_fugacity_mixture(
                model_name=model_name,
                model_input=model_input,
                model_source=model_source
            )

            # NOTE: extract fugacity coefficient for each component
            res_1 = res['vapor']

            # res
            phi_comp = {}

            # looping through each component
            for i, key in enumerate(res_1.keys()):
                # set
                phi_comp[key] = res_1[key]['fugacity_coefficient']['value']

            return phi_comp
        except Exception as e:
            raise Exception(
                f"Error in calculating the fugacity coefficient for the gaseous mixture: {str(e)}") from e

    def _cal_activity_coefficient_solution(
        self,
        model_name: str,
        model_input: Dict[str, Any],
    ):
        """
        Calculate the activity coefficient of solution using the specified model.

        Parameters
        ----------
        model_name : str
            The model to use for calculation. Options: "NRTL", "UNIQUAC".
        model_input : dict
            The input data for the model.
            - mole_fraction: Dictionary of component mole fractions.
            - tau_ij: Dictionary of interaction parameters.
            - alpha_ij: Dictionary of interaction parameters.
            - r_i: relative van der Waals volume of component i
            - q_i: relative surface area of component i
        model_source : dict
            The source of the model data.
            - datasource: Data source for the model.
            - equationsource: Equation source for the model.

        Returns
        -------
        dict
            A dictionary containing the fugacity coefficients for each component in the liquid mixture.
        """
        try:
            # SECTION: model source
            # model_source = {
            #     "datasource": self.datasource,
            #     "equationsource": self.equationsource
            # }
            # check activity
            if self.activity is None or self.activity == 'None':
                raise ValueError(
                    f"Invalid activity model. Must be {ACTIVITY_MODELS}.")

            # NOTE: model name
            if model_name == "NRTL":
                # calculate fugacity
                res, others = self.activity.cal(model_input=model_input)
            elif model_name == "UNIQUAC":
                # calculate fugacity
                res, others = self.activity.cal(model_input=model_input)
            else:
                raise ValueError(
                    f"Invalid model name: {model_name}. Must be 'NRTL' or 'UNIQUAC'.")

            # NOTE: extract fugacity coefficient for each component
            return others['AcCo_i_comp']
        except Exception as e:
            raise Exception(
                f"Error in calculating the fugacity coefficient for the liquid mixture: {str(e)}") from e

    def calc_chemical_potential_pure(
        self,
        temperature: float
    ) -> Dict[str, Any]:
        '''
        Calculate the pure chemical potential of a reaction system at a given temperature.
        The chemical potential of pure component contains the Gibbs energy of formation at 298.15 K and deviation from the standard state at the given temperature. No correction for pressure is applied, and also no correction for non-ideality is applied.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K] at which to calculate the chemical potential.


        Returns
        -------
        Dict[str, Any]
            A dictionary containing the chemical potential for each component at the specified temperature.
        '''
        try:
            # NOTE:
            # SECTION: calculate the chemical potential at the given temperature
            ChePot_PURE = {}

            # loop through each component
            for component in self.component_state_list:
                # NOTE: check if component is valid
                if len(component) != 3:
                    raise ValueError(
                        f"Invalid component state: {component}. It should be a tuple of (molecule, state, molecule_state).")

                # get component name
                component_name = [component[0], component[2]]
                key_ = f"{component[2]}"

                # calculate chemical potential at the given temperature
                res_ = self.ReactionAnalyzer_.component_energy_at_temperature(
                    datasource=self.datasource,
                    equationsource=self.equationsource,
                    component_names=component_name,
                    temperature=temperature,
                    res_format='symbolic',
                )

                # set
                ChePot_PURE[key_] = res_['GiEnFo_T']

            # NOTE: return chemical potential pure at T
            return ChePot_PURE
        except Exception as e:
            raise Exception(
                f"Error in calculating the gibbs energy of formation: {str(e)}") from e

    def chemical_potential_gas_mixture(self,
                                       temperature: Temperature,
                                       pressure: Pressure,
                                       chemical_potential_pure: Dict[str, Any],
                                       gas_mixture: str = 'ideal',
                                       **kwargs
                                       ) -> Dict[str, float]:
        '''
        Calculate the component fugacity coefficient for gaseous mixture.

        Parameters
        ----------
        temperature : Temperature
            Temperature reference for the reaction system.
            - value: float
            - unit: 'K'
        pressure : Pressure
            Pressure reference for the reaction system.
            - value: float
            - unit: 'bar'
        chemical_potential_std : dict
            Dictionary containing the standard chemical potential for each component.
            - key: component name (e.g., 'CO2-g')
            - value: standard chemical potential value (float)
        gas_mixture : str, optional
            The gas mixture model to use for the calculation. Options: 'ideal', 'SRK', 'PR', 'RK'.
        **kwargs : dict, optional
            Additional keyword arguments for the calculation, such as temperature and pressure.

        Returns
        -------
        Dict[str, float] | None
            A dictionary containing the chemical potential for each component in the gaseous mixture.
        '''
        try:
            # SECTION: prepare model input
            # NOTE: gas mixture
            component_list: List[str] = self.phase_contents.get('g', [])
            # NOTE: phase stream
            phase_stream_gas = self.phase_stream.get('g', {})

            # NOTE: continue calculation
            # ! check if component_list is empty
            if len(component_list) == 0:
                return {}

            # NOTE: get mole fraction
            N0s = {}

            # looping through phase stream
            for key, value in phase_stream_gas.items():
                # check if key is in component list
                if key in component_list:
                    # set mole fraction
                    N0s[key] = value['phase_mole_fraction']

            # SECTION: model input
            # set
            model_input = {
                'feed-specification': N0s,
                'pressure': [
                    pressure['value'],
                    pressure['unit']
                ],
                'temperature': [
                    temperature['value'],
                    temperature['unit']
                ],
            }

            # SECTION: calculate the chemical potential
            # init
            ChePot_comp = {}

            # check
            if gas_mixture == 'ideal':
                # fugacity coefficient for ideal gas mixture
                phi_comp = {
                    component: 1.0 for component in component_list
                }

            elif gas_mixture == 'non-ideal':
                # calculate the fugacity coefficient
                phi_comp = self._cal_fugacity_coefficient_gaseous_mixture(
                    model_name=self.eos_model,
                    model_input=model_input
                )
            else:
                raise ValueError(
                    f"Invalid gas mixture model: {gas_mixture}. Must be 'ideal' or 'non-ideal'.")

            # NOTE: calculate the chemical potential for each component
            # looping through each component
            for component in component_list:
                # NOTE: get chemical potential [J/mol]
                ChePot_std_ = chemical_potential_pure.get(component, None)

                # check
                if ChePot_std_ is None:
                    raise ValueError(
                        f"Chemical potential for component {component} not found in standard chemical potential dictionary.")

                # NOTE: calculate chemical potential
                # mole fraction
                y_i = N0s[component]
                phi_i = phi_comp[component]
                T = temperature['value']  # temperature [K]
                P = pressure['value']  # pressure [bar]
                P0 = self.P_Ref_bar  # reference pressure [bar]
                R = self.R  # universal gas constant [J/mol.K]
                # mixture term [J/mol]
                _term_ = y_i * phi_i
                log_term_ = log(_term_)
                # calc
                mixture_term_ = R * T * log_term_

                # calc [J/mol]
                ChePot_comp[component] = ChePot_std_['value'] + mixture_term_

            # NOTE: res
            return ChePot_comp

        except Exception as e:
            raise Exception(
                f"Error in calculating the gas mixture term: {str(e)}") from e

    def chemical_potential_solution(self,
                                    temperature: Temperature,
                                    pressure: Pressure,
                                    chemical_potential_pure: Dict[str, Any],
                                    solution: str = 'ideal',
                                    **kwargs
                                    ) -> Dict[str, float]:
        '''
        Calculate the component chemical potential for solution.

        Parameters
        ----------
        temperature : Temperature
            Temperature reference for the reaction system.
            - value: float
            - unit: 'K'
        pressure : Pressure
            Pressure reference for the reaction system.
            - value: float
            - unit: 'bar'
        chemical_potential_std : dict
            Dictionary containing the standard chemical potential for each component.
            - key: component name (e.g., 'CO2-g')
            - value: standard chemical potential value (float)
        solution : str, optional
            The solution model to use for the calculation. Options: 'ideal', 'NRTL', 'UNIQUAC'.
        **kwargs : dict, optional
            Additional keyword arguments for the calculation, such as temperature and pressure.

        Returns
        '''
        try:
            # SECTION: prepare model input
            # NOTE: gas mixture
            component_list: List[str] = self.phase_contents.get('l', [])
            # NOTE: phase stream
            phase_stream_gas = self.phase_stream.get('l', {})

            # NOTE: continue calculation
            # ! check if component_list is empty
            if len(component_list) == 0:
                return {}

            # NOTE: get mole fraction
            N0s = {}

            # looping through phase stream
            for key, value in phase_stream_gas.items():
                # check if key is in component list
                if key in component_list:
                    # set mole fraction
                    N0s[key] = value['phase_mole_fraction']

            # SECTION: model input
            # set
            model_input = {
                'mole_fraction': N0s,
                'temperature': [
                    temperature['value'],
                    temperature['unit']
                ],
            }

            # SECTION: calculate the chemical potential
            # init
            ChePot_comp = {}

            # check
            if solution == 'ideal':
                # activity coefficient for ideal solution
                AcCo_i_comp = {
                    component: 1.0 for component in component_list
                }

            elif solution == 'non-ideal':
                # calculate the activity coefficient
                AcCo_i_comp = self._cal_activity_coefficient_solution(
                    model_name=self.activity_model,
                    model_input=model_input,
                )
            else:
                raise ValueError(
                    f"Invalid solution model: {solution}. Must be 'ideal', 'NRTL', or 'UNIQUAC'.")

            # NOTE: calculate the chemical potential for each component
            # looping through each component
            for component in component_list:
                # NOTE: get chemical potential [J/mol]
                ChePot_std_ = chemical_potential_pure.get(component, None)

                # check
                if ChePot_std_ is None:
                    raise ValueError(
                        f"Chemical potential for component {component} not found in standard chemical potential dictionary.")

                # NOTE: calculate chemical potential
                # mole fraction
                x_i = N0s[component]
                AcCo_i = AcCo_i_comp[component]
                T = temperature['value']  # temperature [K]
                P = pressure['value']  # pressure [bar]
                R = self.R  # universal gas constant [J/mol.K]
                # mixture term [J/mol]
                _term_ = x_i * AcCo_i
                log_term_ = log(_term_)
                # calc
                mixture_term_ = R * T * log_term_

                # calc [J/mol]
                ChePot_comp[component] = ChePot_std_['value'] + mixture_term_

            # NOTE: res
            return ChePot_comp

        except Exception as e:
            raise Exception(
                f"Error in calculating the solution term: {str(e)}") from e

    def cal_actual_gibbs_energy_of_reaction(
        self,
        operating_conditions: OperatingConditions,
        gas_mixture: str,
        solution: str
    ):
        '''
        Calculate the actual Gibbs energy of reaction at a given temperature and pressure.

        Parameters
        ----------
        temperature : Temperature
            Temperature reference for the reaction system.
            - value: float
            - unit: 'K'
        pressure : Pressure
            Pressure reference for the reaction system.
            - value: float
            - unit: 'bar'
        gas_mixture : str
            Gas mixture model to use for the calculation. Options: 'ideal', 'SRK', 'PR', 'RK'.
            - ideal: Ideal gas mixture model.
            - non-ideal: Non-ideal gas mixture model using SRK, PR, or RK EOS.
        solution : str
            Solution model to use for the calculation. Options: 'ideal', 'NRTL', 'UNIQUAC'.
            - ideal: Ideal solution model.
            - non-ideal: Non-ideal solution model using NRTL or UNIQUAC activity models.

        '''
        try:
            # SECTION: set operating conditions
            # temperature [K]
            temperature = operating_conditions['temperature']
            # pressure [bar]
            pressure = operating_conditions['pressure']

            # SECTION: calculate the standard chemical potential at the given temperature
            chemical_potential_pure_comp = self.calc_chemical_potential_pure(
                temperature=temperature['value'],
            )

            # SECTION: gas chemical potential
            ChePot_MIXTURE = self.chemical_potential_gas_mixture(
                temperature=temperature,
                pressure=pressure,
                chemical_potential_pure=chemical_potential_pure_comp,
                gas_mixture=gas_mixture,
            )

            # SECTION: solution chemical potential
            ChePot_SOLUTION = self.chemical_potential_solution(
                temperature=temperature,
                pressure=pressure,
                chemical_potential_pure=chemical_potential_pure_comp,
                solution=solution,
            )

            # SECTION: combine chemical potentials
            ChePot = {
                **ChePot_MIXTURE,
                **ChePot_SOLUTION
            }

            # NOTE: calculate the actual Gibbs energy for each reaction
            res = {}

            for reaction_name, reaction in self.reaction_analysis.items():
                # calculate the actual Gibbs energy of reaction
                res[reaction_name] = \
                    self.ReactionAnalyzer_\
                    .calc_actual_gibbs_energy_of_reaction(
                        reaction=reaction,
                        chemical_potential=ChePot,
                        temperature=temperature,
                        pressure=pressure,
                )

                # res
            return {}
        except Exception as e:
            raise Exception(
                f"Error in calculating the chemical potential at the given temperature: {str(e)}") from e
