# import libs
from typing import List, Dict, Any, Literal, Optional
import pycuc
# local
from .reactionanalyzer import ReactionAnalyzer
from ..utils import ChemReactUtils


class Reaction():
    """Class to represent a chemical reaction."""

    # variables
    reaction_analysis_result = {}

    def __init__(self, datasource, equationsource, reaction):
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
        '''
        # set
        self.datasource = datasource
        self.equationsource = equationsource
        self.reaction = reaction

        # NOTE: init
        self.ChemReactUtils_ = ChemReactUtils()
        self.ReactionAnalyzer_ = ReactionAnalyzer()

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
            _res_0 = self.ChemReactUtils_.analyze_reaction(self.reaction)

            # NOTE: energy analysis
            # set input
            _input = {**reaction_res, **_res_0}

            energy_analysis_result = self.ReactionAnalyzer_.energy_analysis(
                self.datasource,
                self.equationsource,
                _input)

            # update
            self.name = name
            self.reaction = reaction
            self.reactants = _res_0['reactants']
            self.products = _res_0['products']
            self.reaction_coefficient = _res_0['reaction_coefficient']
            self.carbon_count = _res_0['carbon_count']
            self.energy_analysis_result = energy_analysis_result
            # extract
            self.GiEnFo = energy_analysis_result['GiEnFo']
            self.EnFo = energy_analysis_result['EnFo']
            self.GiEn_rxn_298 = energy_analysis_result['GiEn_rxn_298']
            self.En_rxn_298 = energy_analysis_result['En_rxn_298']
            self.K_eq = energy_analysis_result['K_eq']

            # reaction analysis result
            self.reaction_analysis_result = {
                'name': name,
                'reaction': reaction,
                'reactants': _res_0['reactants'],
                'products': _res_0['products'],
                'reaction_coefficient': _res_0['reaction_coefficient'],
                'carbon_count': _res_0['carbon_count'],
                'energy_analysis_result': energy_analysis_result
            }

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.go(): {str(e)}") from e

    def Keq_T(self, temperature: List[float | str]) -> Dict[str, Any]:
        """
        Calculate the equilibrium constant at a given temperature.

        Parameters
        ----------
        temperature : list[float, str]
            Temperature in any unit, e.g. [300.0, "K"].
            It is automatically converted to Kelvin.

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

            # NOTE: calculate equilibrium constant at a given temperature
            res = self.ReactionAnalyzer_.vh(
                self.datasource,
                self.equationsource,
                T,
                self.reaction_analysis_result,
            )

            return res
        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.Keq_T(): {str(e)}") from e
