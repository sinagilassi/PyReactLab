# import libs
from typing import Dict, Any, List, Literal, Optional
# local
from .reaction import Reaction
from .thermolinkdb import ThermoLinkDB
from .refmanager import ReferenceManager
from ..configs import DATASOURCE, EQUATIONSOURCE


class ReactionSystem(ThermoLinkDB):
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
        self.system_source = model_source

        self.__model_source = model_source
        self.__datasource = {
        } if model_source is None else model_source[DATASOURCE]
        self.__equationsource = {
        } if model_source is None else model_source[EQUATIONSOURCE]

        # NOTE: init class
        ThermoLinkDB.__init__(self)
        ReferenceManager.__init__(self)

        # NOTE: load reference
        # reference plugin
        self._references = self.load_reference()

        # NOTE: component datasource and equationsource
        self._datasource, self._equationsource = self.link_thermodb()

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

    def link_thermodb(self) -> None:
        """
        Link the reaction system to a thermodynamic database.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        try:
            # SECTION: set datasource and equationsource
            # NOTE: check if datasource and equationsource are provided in system_source
            # datasource
            datasource = self.system_source.get('datasource')
            # equationsource
            equationsource = self.system_source.get('equationsource')

            # check if datasource and equationsource are provided
            if datasource is None or equationsource is None:
                raise ValueError(
                    "Datasource and equationsource must be provided in system_source.")

            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # SECTION:
            # build datasource
            component_datasource = self.set_datasource(
                self.components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                self.components, reference)

        except Exception as e:
            raise Exception(
                f"Error linking thermodynamic database: {e}") from e
