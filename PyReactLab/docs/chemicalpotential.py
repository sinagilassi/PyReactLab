# import libs
import numpy as np
from typing import Literal, Dict, Any, List
# local
# local
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K, EOS_MODELS, ACTIVITY_MODELS
)
from .reaction import Reaction


class ChemicalPotential:
    """
    Class to calculate chemical potential of a reaction system.
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
                 comp_list: List[Dict[str, float | int]],
                 reaction_analysis: Dict,
                 phase_stream: Dict[str, Any],
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
        # reaction list
        self.reaction_list = reaction_list
        # set component list
        self.comp_list = comp_list
        # set reaction analysis
        self.reaction_analysis = reaction_analysis
        # phase stream
        self.phase_stream = phase_stream
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
        return self._eos_model

    @eos_model.setter
    def eos_model(self, value: str):
        '''Set eos model'''
        # set
        self._eos_model = value

    @property
    def activity_model(self):
        '''Get activity model name'''
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
        model_name: Literal[
            "SRK", "PR", "RK"
        ],
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
        model_name: Literal[
            'NRTL', 'UNIQUAC'
        ],
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

    def cal_actual_gibbs_energy_of_reaction(
        self,
        phase_contents: Dict[str, Any],
        phase_stream: Dict[str, Any],
        temperature: float,
        pressure: float,
        gas_mixture: str = 'ideal',
        solution: str = 'ideal',
    ):
        pass
