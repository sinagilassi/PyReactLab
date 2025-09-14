# import libs
import logging
from typing import List, Dict, Optional, Any
import pycuc
from pyThermoDB.core import (
    TableData,
    TableEquation,
    TableMatrixData
)
# local
from ..configs import DATASOURCE, EQUATIONSOURCE

# NOTE: get logger
logger = logging.getLogger(__name__)


class ThermoSource:
    '''
    Source class for the pyThermoFlash package.
    '''
    # NOTE: variables

    def __init__(
        self,
        model_source: Optional[Dict] = None,
        **kwargs
    ):
        '''Initialize the Source class.'''
        # set
        self.model_source = model_source

        # NOTE: source
        if model_source is None:
            self._datasource = None
            self._equationsource = None
        else:
            # set
            self._datasource, self._equationsource = self.set_source(
                model_source=model_source
            )

    def __repr__(self):
        des = "Source class for the pyThermoFlash package."
        return des

    @property
    def datasource(self) -> Dict:
        '''
        Get the datasource property.

        Returns
        -------
        dict
            The datasource dictionary.
        '''
        # NOTE: check if model source is valid
        if self._datasource is None:
            return {}
        return self._datasource

    @property
    def equationsource(self) -> Dict:
        '''
        Get the equationsource property.

        Returns
        -------
        dict
            The equationsource dictionary.
        '''
        # NOTE: check if model source is valid
        if self._equationsource is None:
            return {}
        return self._equationsource

    def set_source(self, model_source: Dict):
        '''
        Set the model source.

        Parameters
        ----------
        model_source : dict
            The model source dictionary.
        '''
        # NOTE: source
        # datasource
        _datasource = {} \
            if model_source is None else model_source[DATASOURCE]

        # equationsource
        _equationsource = {}\
            if model_source is None else model_source[EQUATIONSOURCE]

        # res
        return _datasource, _equationsource

    def eq_extractor(
        self,
        component_name: str,
        prop_name: str
    ) -> Dict[str, TableEquation]:
        '''
        Extracts the equilibrium property from the model source.

        Parameters
        ----------
        component_name : str
            The name of the component.
        prop_name : str
            The name of the property to extract.

        Returns
        -------
        dict[str, TableEquation]
            The extracted property.
        '''
        try:
            if self.equationsource is None:
                logger.warning("Equation source is not defined.")
                return {}

            # NOTE: check component
            if component_name not in self.equationsource.keys():
                logger.error(
                    f"Component '{component_name}' not found in model source.")
                return {}

            # NOTE: check property
            if prop_name not in self.equationsource[component_name].keys():
                logger.error(
                    f"Property '{prop_name}' not found in model source registered for {component_name}.")
                return {}

            # return
            return self.equationsource[component_name][prop_name]
        except Exception as e:
            logger.error(f"Error extracting equation: {e}")
            raise

    def data_extractor(
        self,
        component_name: str,
        prop_name: str
    ) -> Dict:
        '''
        Extracts the data property from the model source.

        Parameters
        ----------
        component_name : str
            The name of the component.
        prop_name : str
            The name of the property to extract.

        Returns
        -------
        dict
            The extracted property.

        Notes
        -----
        The extracted property is a dictionary containing the following keys:
        - 'symbol': The symbol of the property.
        - 'property_name': The name of the property.
        - 'unit': The unit of the property.
        - 'value': The value of the property.
        - 'message': A message about the property.
        '''
        try:
            if self.datasource is None:
                logger.warning("Data source is not defined.")
                return {}

            # NOTE: check component
            if component_name not in self.datasource.keys():
                logger.error(
                    f"Component '{component_name}' not found in model datasource.")
                return {}

            # NOTE: check property
            if prop_name not in self.datasource[component_name].keys():
                logger.error(
                    f"Property '{prop_name}' not found in model datasource registered for {component_name}.")
                return {}

            # return
            return self.datasource[component_name][prop_name]
        except Exception as e:
            logger.error(f"Error extracting data: {e}")
            raise

    def check_args(
        self,
        component_name: str,
        args
    ):
        '''
        Checks equation args

        Parameters
        ----------
        component_name : str
            The name of the component.
        args : tuple
            equation args
        '''
        try:
            # required args
            required_args = []

            # datasource list
            datasource_component_list = list(
                self.datasource[component_name].keys())

            # NOTE: default args
            datasource_component_list.append("P")
            datasource_component_list.append("T")

            # check args within datasource
            for arg_key, arg_value in args.items():
                # symbol
                if arg_value['symbol'] in datasource_component_list:
                    # update
                    required_args.append(arg_value)
                else:
                    raise Exception('Args not in datasource!')

            # res
            return required_args

        except Exception as e:
            raise Exception('Finding args failed!, ', e)

    def build_args(
        self,
        component_name: str,
        args,
        ignore_symbols: List[str] = ["T", "P"]
    ):
        '''
        Builds args

        Parameters
        ----------
        component_name : str
            The name of the component.
        args : tuple
            equation args
        ignore_symbols : list
            list of symbols to ignore, default is ["T", "P"]
        '''
        try:
            # res
            res = {}
            for arg in args:
                # symbol
                symbol = arg['symbol']

                # NOTE: check if symbol is in ignore symbols
                if ignore_symbols is not None:
                    # check in ignore symbols
                    if symbol not in ignore_symbols:
                        # check in component database
                        for key, value in self.datasource.items():
                            if symbol == key:
                                res[symbol] = value
                else:
                    # check in component database
                    for key, value in self.datasource.items():
                        if symbol == key:
                            res[symbol] = value
            return res
        except Exception as e:
            raise Exception('Building args failed!, ', e)

    def eq_builder(
        self,
        components: List[str],
        prop_name: str,
        **kwargs
    ):
        '''
        Builds the equation for the given components and property name.

        Parameters
        ----------
        components : List[str]
            List of component names to build the equation for.
        prop_name : str
            The name of the property to build the equation for.
        **kwargs : dict
            Additional keyword arguments for the equation builder.

        Returns
        -------
        dict
            The built equation for each component.
        '''
        # NOTE: check if model source is valid
        if self.equationsource is None:
            raise ValueError("Equation source is not defined.")

        # NOTE: check property
        for component in components:
            # check equation availability
            if prop_name not in self.equationsource[component].keys():
                raise ValueError(
                    f"Property '{prop_name}' not found in model source registered for {component}.")

        # NOTE: property name
        if len(prop_name) == 0:
            raise ValueError("Property name cannot be empty.")

        # strip
        prop_name = prop_name.strip()

        # NOTE: vapor pressure source
        eq_src_comp = {}

        # looping through components
        for component in components:
            # NOTE: equation source
            _eq = None
            # select equation [?]
            _eq = self.eq_extractor(component, prop_name)

            # check
            if _eq is None:
                raise ValueError(
                    f"Equation for property '{prop_name}' not found for component '{component}'.")

            # check type
            if not isinstance(_eq, TableEquation):
                raise ValueError(
                    f"Equation for property '{prop_name}' for component '{component}' is not a valid TableEquation instance.")

            # NOTE: args
            _args = _eq.args
            # check args (SI)
            _args_required = self.check_args(
                component, _args)

            # build args
            _args_ = self.build_args(
                component, _args_required)

            # NOTE: update P and T
            _args_['T'] = None
            _args_['P'] = None

            # set
            eq_src_comp[component] = {
                "value": _eq,
                "args": _args_,
                "return": _eq.returns
            }

        # res
        return eq_src_comp

    def exec_eq(
        self,
        components: List[str],
        eq_src_comp: Dict[str, Any],
        prop_res_unit: str,
        args_values: Optional[Dict[str, float]] = None,
        **kwargs
    ):
        '''

        Executes the equation for the given components and arguments.

        Parameters
        ----------
        components : List[str]
            List of component names to execute the equation for.
        eq_src_comp : Dict[str, Any]
            Dictionary containing the equation source for each component.
        args_values : Dict[str, float]
            Dictionary containing the values for the arguments.
        **kwargs : dict
            Additional keyword arguments for the equation execution.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing the results of the equation execution for each component.
        '''
        try:
            # NOTE: check if model source is valid
            if not isinstance(components, list):
                raise ValueError("Components must be a list.")

            # NOTE: check if eq_src_comp is a dictionary
            if not isinstance(eq_src_comp, dict):
                raise ValueError("Equation source must be a dictionary.")

            # NOTE: check if args_values is a dictionary
            if args_values is not None and not isinstance(args_values, dict):
                raise ValueError("Arguments values must be a dictionary.")

            # NOTE: check if components are in eq_src_comp
            for component in components:
                if component not in eq_src_comp:
                    raise ValueError(
                        f"Component '{component}' not found in equation source.")

            # NOTE: check if prop_res_unit is a string
            if not isinstance(prop_res_unit, str):
                raise ValueError("Property result unit must be a string.")

            # strip
            prop_res_unit = prop_res_unit.strip()

            # NOTE: init res
            prop_res = []
            prop_res_dict = {}

            # looping over components
            for i, component in enumerate(components):
                # NOTE: equation [unit:?]
                eq_ = eq_src_comp[component]['value']
                args_ = eq_src_comp[component]['args']

                # update args
                # args_['T'] = T
                if args_values is not None:
                    # NOTE: check if args_values is a dictionary
                    if not isinstance(args_values, dict):
                        raise ValueError(
                            "Arguments values must be a dictionary.")

                    # NOTE: update args with args_values
                    for key, value in args_values.items():
                        # NOTE: check if key is in args
                        if key in args_:
                            # update
                            args_[key] = value

                # NOTE: cal
                res_ = eq_.cal(**args_)

                # extract
                res_value_ = res_['value']
                res_unit_ = res_['unit']

                # NOTE: convert to Pa
                unit_block_ = f"{res_unit_} => {prop_res_unit}"
                prop_value = pycuc.to(res_value_, unit_block_)

                # NOTE: return
                return_key_ = []
                return_val_ = []
                for key, value in eq_.returns.items():
                    return_key_.append(key)
                    return_val_.append(value)

                # check length
                property_name = ''
                for key in return_key_:
                    if len(property_name) > 0:
                        property_name += ', '
                    property_name += key

                property_symbol = ''
                for val in return_val_:
                    if len(property_symbol) > 0:
                        property_symbol += ', '
                    property_symbol += str(val['symbol'])

                # save
                prop_res.append(prop_value)
                # dict type
                prop_res_dict[component] = {
                    "value": prop_value,
                    "unit": prop_res_unit,
                    "symbol": property_symbol,
                    "property_name": property_name
                }

            # NOTE: results
            return prop_res, prop_res_dict
        except Exception as e:
            raise Exception('Executing equation failed!, ', e)

    def get_component_data(
            self,
            component_name: str,
            components: List[str]
    ):
        """
        Get the component data from the datasource.

        Parameters
        ----------
        component_name : str
            The name of the component.

        Returns
        -------
        dict
            The component data.
        """
        try:
            # check
            if not isinstance(component_name, str):
                raise ValueError("Component name must be a string.")

            # check available
            if not isinstance(components, list):
                raise ValueError("Components must be a list.")

            if component_name not in components:
                raise ValueError(
                    f"Component {component_name} is not available in the system.")

            # data
            data = {}

            # check datasource
            if (self.datasource is not None and
                    isinstance(self.datasource, dict)):
                # add
                dt_ = self.datasource.get(component_name)
                if dt_ is not None:
                    data.update(dt_)

            # check equationsource
            if (self.equationsource is not None and
                    isinstance(self.equationsource, dict)):
                # add
                eq_ = self.equationsource.get(component_name)
                if eq_ is not None:
                    data.update(eq_)

            # get data
            return data
        except Exception as e:
            raise Exception(f"Error in component_data: {e}") from e
