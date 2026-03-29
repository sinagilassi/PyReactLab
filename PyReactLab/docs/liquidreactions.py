# Liquid Reactions
# import libs
import logging
from typing import List, Dict, Any
from pyThermoLinkDB.models import ModelSource
from pyThermoLinkDB.thermo import Source
# local
from ..core.reactionsystem import ReactionSystem

# NOTE: logger
logger = logging.getLogger(__name__)


class LiquidReactionSystem(ReactionSystem):
    # NOTE: variables
    # phase
    _phase: str = 'liquid'

    def __init__(
        self,
        system_name: str,
        reactions: List[Dict[str, Any]],
        model_source: ModelSource,
        source: Source,
    ):
        """
        Initialize the LiquidReactionSystem class.

        Parameters
        ----------
        system_name : str
            Name of the reaction system.
        reactions : list
            List of reactions in the system must be in the form of a list of dictionaries as the following keys
            - 'reaction': str, the reaction equation.
            - 'name': str, the name of the reaction.
        model_source : ModelSource
            Inputs for the reaction system which contains the data source and equation source for the components in the reaction system.
        source : Source
            The source object containing the model source and component key for the reaction system.
        """
        super().__init__(
            system_name=system_name,
            reactions=reactions,
            model_source=model_source,
            source=source,
            phase_rule=self._phase
        )

        # NOTE: check phase
        self.checking_phase()

    @property
    def phase(self) -> str:
        """
        Get the phase of the reaction system.

        Returns
        -------
        str
            Phase of the reaction system.
        """
        return self._phase

    def checking_phase(self) -> None:
        """
        Check the phase of the reaction system.
        """
        try:
            # overall reaction phase
            overall_reaction_phase = self.overall_reaction_phase
            # check overall reaction phase
            if overall_reaction_phase:
                if isinstance(overall_reaction_phase, str):
                    if overall_reaction_phase.lower() != self._phase:
                        raise ValueError(
                            f"Overall reaction phase is {overall_reaction_phase} but expected {self._phase}.")
                else:
                    raise ValueError(
                        f"Overall reaction phase is {overall_reaction_phase} but expected {self._phase}.")
            else:
                raise ValueError(
                    f"Overall reaction phase is {overall_reaction_phase} but expected {self._phase}.")
        except Exception as e:
            raise Exception(f"Error checking phase: {e}") from e
