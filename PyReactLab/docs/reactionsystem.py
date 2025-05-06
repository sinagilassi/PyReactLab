# import libs
from typing import Dict, Any, List, Literal, Optional
# local
from .reaction import Reaction
from .thermolinkdb import ThermoLinkDB
from .refmanager import ReferenceManager
from .reactionanalyzer import ReactionAnalyzer
from ..utils import ChemReactUtils


class ReactionSystem(ThermoLinkDB, ReferenceManager):
    """Class to represent a system of chemical reactions."""

    # NOTE: class variables
    __system_name = None
    __reactions = None

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

        # SECTION:

    def __str__(self):
        """String representation of the reaction system."""
        return "\n".join([str(reaction) for reaction in self.__reactions])

    @property
    def system_name(self) -> str:
        """Get the name of the reaction system."""
        return self.__system_name

    @property
    def reactions(self) -> List[Dict[str, Any]]:
        """Get the reactions of the reaction system."""
        return self.__reactions

    def go(self) -> None:
        """
        Execute the primary analysis for the reaction system.
        """
        # NOTE: initialize
        ChemReactUtils_ = ChemReactUtils()
        ReactionAnalyzer_ = ReactionAnalyzer()

        # SECTION: reaction system analysis
        # analyze reaction
        reaction_res = {}
        for item in self.reactions:
            _res = ChemReactUtils_.analyze_reaction(item)
            # name
            name = item['name']
            reaction_res[name] = _res

        # SECTION: analyze overall reaction
        res_0 = ChemReactUtils_.analyze_overall_reactions(
            self.reactions)

        # SECTION: set component
        res_1 = ChemReactUtils_.define_component_id(
            self.reaction_res)
        # extract
        component_list, component_dict, comp_list, comp_coeff = res_1

        # SECTION: energy analysis
        # energy analysis result list
        self.energy_analysis_res_list = {}

        # loop through each reaction
        for item in reaction_res:
            _res = ReactionAnalyzer_.energy_analysis(
                self.datasource,
                self.equationsource,
                reaction_res[item])

            # save
            self.energy_analysis_res_list[item] = _res
