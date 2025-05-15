# Gas Reactions
# import libs
from typing import List, Dict, Any
# local
from .reactionsystem import ReactionSystem


class GasReactionSystem(ReactionSystem):
    # NOTE: variables
    # phase
    _phase: str = 'gas'

    def __init__(self,
                 system_name: str,
                 reactions: List[Dict[str, Any]],
                 model_source: Dict[str, Any]
                 ):
        """
        Initialize the GasReactionSystem class.

        Parameters
        ----------
        system_name : str
            Name of the reaction system.
        reactions : list
            List of reactions in the system must be in the form of a list of dictionaries as the following keys
            - 'reaction': str, the reaction equation.
            - 'name': str, the name of the reaction.
        model_source : dict
            Inputs for the reaction system which
        """
        super().__init__(system_name=system_name,
                         reactions=reactions,
                         model_source=model_source,
                         phase_rule=self._phase)

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
