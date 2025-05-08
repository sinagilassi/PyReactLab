# import libs
from typing import Dict, Any, List, Literal, Optional
import pycuc
# local
from .reaction import Reaction
from .thermolinkdb import ThermoLinkDB
from .refmanager import ReferenceManager
from .reactionanalyzer import ReactionAnalyzer
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
                # _res = ChemReactUtils_.analyze_reaction(item)

                # name
                name = item['name']
                # update
                reaction_res[name] = _res
                # set
                self.__reaction_list[name] = r_

            # SECTION: analyze overall reaction
            res_0 = ChemReactUtils_.analyze_overall_reactions(
                self.reactions)

            # SECTION: set component
            res_1 = ChemReactUtils_.define_component_id(
                reaction_res)
            # extract
            component_list, component_dict, comp_list, comp_coeff = res_1

            # SECTION: energy analysis
            # energy analysis result list
            energy_analysis_res_list = {}

            # loop through each reaction
            for item in reaction_res:
                # name
                name = item

                # NOTE: energy analysis
                # _res = ReactionAnalyzer_.energy_analysis(
                #     self.datasource,
                #     self.equationsource,
                #     reaction_res[item])

                reaction_ = self.__reaction_list[name]
                _res = reaction_.energy_analysis_result

                # save
                energy_analysis_res_list[item] = _res

            # NOTE: set primary analysis result
            self.reaction_res = reaction_res
            self.overall_reaction_res = res_0
            self.component_list = component_list
            self.component_dict = component_dict
            self.comp_list = comp_list
            self.comp_coeff = comp_coeff
            self.energy_analysis_res_list = energy_analysis_res_list

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.go(): {str(e)}") from e

    def equilibrium_constant_at_temperature(self,
                                            reaction_name: str,
                                            temperature: list[float | str]):
        """
        Calculate the equilibrium constant at a given temperature.

        Parameters
        ----------
        reaction_name : str
            Name of the reaction.
        temperature : float
            Temperature in Kelvin.

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
            return reaction.Keq_T(temperature)

        except Exception as e:
            raise Exception(
                f"Error in ReactionSystem.equilibrium_constant_at_temperature(): {str(e)}") from e
