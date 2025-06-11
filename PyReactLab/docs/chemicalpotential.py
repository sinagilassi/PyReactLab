# import libs
from typing import Literal
# local
from .reactionanalyzer import ReactionAnalyzer
from ..utils import ChemReactUtils


class ChemicalPotential:
    """
    Class to calculate chemical potential of a reaction system.
    """

    def __init__(
        self,
        datasource: dict,
        equationsource: dict,
        reaction,
        **kwargs
    ):
        """
        Initialize the ChemicalPotential class.

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
        """
        # set
        self.datasource = datasource
        self.equationsource = equationsource
        self.reaction = reaction
        self.kwargs = kwargs

        # NOTE: init
        self.ChemReactUtils_ = ChemReactUtils()
        self.ReactionAnalyzer_ = ReactionAnalyzer()

        # NOTE: kwargs
        # phase rule
        self.phase_rule = kwargs.get('phase_rule', None)

        # NOTE: reaction analyzer

    def __mu_gas_pure_ideal(self, mu_std: float):
        pass

    def __mu_gas_pure_non_ideal(self):
        pass

    def __mu_liquid_pure_ideal(self):
        pass

    def __mu_liquid_pure_non_ideal(self):
        pass

    def __mu_gas_mixture_ideal(self):
        pass

    def __mu_gas_mixture_non_ideal(self):
        pass

    def __mu_solution_ideal(self):
        pass

    def __mu_solution_non_ideal(self):
        pass
