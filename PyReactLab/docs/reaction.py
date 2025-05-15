# import libs
from typing import List, Dict, Any, Literal, Optional
import pycuc
# local
from .reactionanalyzer import ReactionAnalyzer
from ..utils import ChemReactUtils
from ..configs import (
    EQUILIBRIUM_CONSTANT_STD, EQUILIBRIUM_CONSTANT_STD_SYMBOL,
    EQUILIBRIUM_CONSTANT, EQUILIBRIUM_CONSTANT_SYMBOL,
    GIBBS_FREE_ENERGY_OF_REACTION_STD, GIBBS_FREE_ENERGY_OF_REACTION_STD_SYMBOL,
    ENTHALPY_OF_REACTION_STD, ENTHALPY_OF_REACTION_STD_SYMBOL,
    GIBBS_FREE_ENERGY_OF_FORMATION_STD, GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL,
    ENTHALPY_OF_FORMATION_STD, ENTHALPY_OF_FORMATION_STD_SYMBOL,
    GIBBS_FREE_ENERGY_OF_REACTION_T, ENTHALPY_OF_REACTION_T,
    GIBBS_FREE_ENERGY_OF_REACTION_T_SYMBOL, ENTHALPY_OF_REACTION_T_SYMBOL,
)


class Reaction():
    """Class to represent a chemical reaction."""
    # NOTE: variables

    # variables
    reaction_analysis_result = {}

    def __init__(self,
                 datasource,
                 equationsource,
                 reaction,
                 **kwargs):
        '''
        Initialize the reaction object.

        Parameters
        ----------
        datasource : dict
            Dictionary containing the thermodynamic data for the reaction.
        equationsource : dict
            Dictionary containing the equations for the reaction.
        reaction : dict
            Reaction object in the form of a dictionary with the following keys:
            - 'name': str, the name of the reaction.
            - 'reaction': str, the reaction equation.
        **kwargs : dict
            Additional keyword arguments.
        '''
        # set
        self.datasource = datasource
        self.equationsource = equationsource
        self.reaction = reaction

        # NOTE: init
        self.ChemReactUtils_ = ChemReactUtils()
        self.ReactionAnalyzer_ = ReactionAnalyzer()

        # NOTE: kwargs
        # phase rule
        self.phase_rule = kwargs.get('phase_rule', None)

        # NOTE: reaction analyzer
        self._reaction_analyzer()

    def _reaction_analyzer(self) -> None:
        """
        Execute the primary analysis for the reaction system.
        """
        try:
            # SECTION: reaction system analysis
            # analyze reaction
            reaction_res = {}

            # NOTE: name
            name = self.reaction['name']
            # update
            reaction_res['name'] = name

            # NOTE: reaction
            reaction = self.reaction['reaction']
            # update
            reaction_res['reaction'] = reaction

            # NOTE: extract reaction information
            # 'name': name,
            # 'reaction': reaction,
            # 'reactants': reactants,
            # 'products': products,
            # 'reaction_coefficient': reaction_coefficient,
            # 'carbon_count': carbon_count
            # 'reaction_state': reaction_state,
            _res_0 = self.ChemReactUtils_.analyze_reaction(
                self.reaction,
                phase_rule=self.phase_rule,
            )

            # NOTE: energy analysis
            # set input
            _input = {**reaction_res, **_res_0}

            energy_analysis = self.ReactionAnalyzer_.energy_analysis(
                self.datasource,
                self.equationsource,
                _input)

            # update
            self.reaction_name = name
            self.reaction_equation = reaction
            self.reactants = _res_0['reactants']
            self.reactants_names = _res_0['reactants_names']
            self.products = _res_0['products']
            self.products_names = _res_0['products_names']
            self.reaction_coefficient = _res_0['reaction_coefficient']
            self.carbon_count = _res_0['carbon_count']
            self.reaction_state = _res_0['reaction_state']
            self.reaction_phase = _res_0['reaction_phase']
            self.state_count = _res_0['state_count']
            self.energy_analysis = energy_analysis
            # extract
            self.gibbs_free_energy_of_formation_std = energy_analysis[
                GIBBS_FREE_ENERGY_OF_FORMATION_STD]
            self.enthalpy_of_formation_std = energy_analysis[
                ENTHALPY_OF_FORMATION_STD
            ]
            self.gibbs_free_energy_of_reaction_std = energy_analysis[
                GIBBS_FREE_ENERGY_OF_REACTION_STD]
            self.enthalpy_of_reaction_std = energy_analysis[
                ENTHALPY_OF_REACTION_STD
            ]
            self.equilibrium_constant_std = energy_analysis[
                EQUILIBRIUM_CONSTANT_STD
            ]

            # reaction analysis result
            self.reaction_analysis_result: Dict[str, Any] = {
                'name': name,
                'reaction': reaction,
                'reactants': _res_0['reactants'],
                'reactants_names': _res_0['reactants_names'],
                'products': _res_0['products'],
                'products_names': _res_0['products_names'],
                'reaction_coefficient': _res_0['reaction_coefficient'],
                'carbon_count': _res_0['carbon_count'],
                'reaction_state': _res_0['reaction_state'],
                'state_count': _res_0['state_count'],
                'energy_analysis': energy_analysis
            }

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.go(): {str(e)}") from e

    def check_state(self):
        """
        Check the state of the reaction system whether it is a liquid, gas, aqueous, or solid or a combination of them.
        """
        try:
            # NOTE: set state
            states = []

            # looping through reactants and products
            for reactant in self.reactants:
                _res = reactant['state']
                # ? add state
                states.append(_res)
            for product in self.products:
                _res = product['state']
                # ? add state
                states.append(_res)

            # NOTE: remove duplicates
            states = list(set(states))
        except Exception as e:
            raise Exception(
                f"Failing in finding the reaction state: {str(e)}") from e

    def cal_equilibrium_constant(self,
                                 temperature: List[float | str],
                                 method: Literal[
                                     "van't Hoff", "shortcut van't Hoff"
                                 ] = "van't Hoff") -> Dict[str, Any]:
        """
        Calculate the equilibrium constant at a given temperature using the van't Hoff equation.

        Parameters
        ----------
        temperature : list[float, str]
            Temperature in any unit, e.g. [300.0, "K"], It is automatically converted to Kelvin.
        method : str, optional
            Method to calculate the equilibrium constant. The default is "van't Hoff".
            Options are "van't Hoff" or "shortcut van't Hoff".

        Returns
        -------
        dict
            Dictionary containing the equilibrium constant at the given temperature.
        """
        try:
            # NOTE: check if T
            # check if T is a list
            if not isinstance(temperature, list):
                raise ValueError("Temperature must be a list.")

            # check if T is a number
            if not isinstance(temperature[0], (int, float)):
                raise ValueError("Temperature must be a number.")

            # check if T is a string
            if not isinstance(temperature[1], str):
                raise ValueError("Temperature unit must be a string.")

            # NOTE: convert temperature to Kelvin
            # set unit
            unit_set = f"{temperature[1]} => K"
            T = pycuc.to(temperature[0], unit_set)

            # SECTION: calculate equilibrium constant at a given temperature
            # check if method is valid
            if method not in ["van't Hoff", "shortcut van't Hoff"]:
                raise ValueError(
                    f"Invalid method: {method}, must be 'van't Hoff' or 'shortcut van't Hoff'.")

            if method == "van't Hoff":
                # ! NOTE: van't Hoff method
                res = self.ReactionAnalyzer_.vh(
                    self.datasource,
                    self.equationsource,
                    T,
                    self.reaction_analysis_result,
                )

                # update
                res['method'] = method
                # return
                return res
            elif method == "shortcut van't Hoff":
                # ! NOTE: shortcut van't Hoff method
                res = self.ReactionAnalyzer_.vh_shortcut(
                    self.datasource,
                    self.equationsource,
                    T,
                    self.enthalpy_of_reaction_std['value'],
                    self.equilibrium_constant_std['value'],
                    self.reaction_name,
                    self.reaction_equation
                )

                # update
                res['method'] = method
                # return
                return res
            else:
                raise ValueError(
                    f"Invalid method: {method}, must be 'van't Hoff' or 'shortcut van't Hoff'.")

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.Keq_T(): {str(e)}") from e

    def cal_reaction_energy(self, temperature: List[float | str]) -> Dict[str, Any]:
        """
        Calculate the reaction energy at a given temperature which consists of the following:
        - Gibbs free energy of reaction at T
        - Enthalpy of reaction at T

        Parameters
        ----------
        temperature : list[float, str]
            Temperature in any unit, e.g. [300.0, "K"], It is automatically converted to Kelvin.

        Returns
        -------
        dict
            Dictionary containing the reaction energy at the given temperature.
        """
        try:
            # NOTE: check if T
            # check if T is a list
            if not isinstance(temperature, list):
                raise ValueError("Temperature must be a list.")

            # check if T is a number
            if not isinstance(temperature[0], (int, float)):
                raise ValueError("Temperature must be a number.")

            # check if T is a string
            if not isinstance(temperature[1], str):
                raise ValueError("Temperature unit must be a string.")

            # NOTE: convert temperature to Kelvin
            # set unit
            unit_set = f"{temperature[1]} => K"
            T = pycuc.to(temperature[0], unit_set)

            # SECTION: calculate reaction energy at a given temperature
            res = self.ReactionAnalyzer_.reaction_energy_analysis(
                self.datasource,
                self.equationsource,
                T,
                self.reaction_analysis_result,
            )

            # NOTE: extract reaction energy
            return {
                "reaction_name": self.reaction_name,
                "reaction_equation": self.reaction_equation,
                "temperature": {
                    "value": T,
                    "symbol": "T",
                    "unit": "K",
                },
                GIBBS_FREE_ENERGY_OF_REACTION_T: res[GIBBS_FREE_ENERGY_OF_REACTION_T],
                ENTHALPY_OF_REACTION_T: res[ENTHALPY_OF_REACTION_T],
            }

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.Keq_T(): {str(e)}") from e
