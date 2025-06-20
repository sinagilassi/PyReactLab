# import libs
from typing import Dict, Any, List, Literal, Optional
from math import exp
from pyThermoDB import TableEquation, TableMatrixEquation, TableData, TableMatrixData
import pycuc
from scipy import integrate
#  local
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K,
    EQUILIBRIUM_CONSTANT_STD, EQUILIBRIUM_CONSTANT_STD_SYMBOL,
    EQUILIBRIUM_CONSTANT, EQUILIBRIUM_CONSTANT_SYMBOL,
    GIBBS_FREE_ENERGY_OF_REACTION_STD, GIBBS_FREE_ENERGY_OF_REACTION_STD_SYMBOL,
    ENTHALPY_OF_REACTION_STD, ENTHALPY_OF_REACTION_STD_SYMBOL,
    GIBBS_FREE_ENERGY_OF_FORMATION_STD, GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL,
    ENTHALPY_OF_FORMATION_STD, ENTHALPY_OF_FORMATION_STD_SYMBOL,
    GIBBS_FREE_ENERGY_OF_REACTION_T, ENTHALPY_OF_REACTION_T,
    GIBBS_FREE_ENERGY_OF_REACTION_T_SYMBOL, ENTHALPY_OF_REACTION_T_SYMBOL,
    GIBBS_FREE_ENERGY_OF_FORMATION_T, ENTHALPY_OF_FORMATION_T,
    GIBBS_FREE_ENERGY_OF_FORMATION_T_SYMBOL, ENTHALPY_OF_FORMATION_T_SYMBOL,
    CHEMICAL_POTENTIAL_MIXTURE_T_P, CHEMICAL_POTENTIAL_MIXTURE_T_P_SYMBOL,
    ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION, ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION_SYMBOL
)
from ..utils import (
    Temperature,
    Pressure
)


class ReactionAnalyzer:
    """Class to analyze a reaction system."""
    # NOTE: class variables
    # universal gas constant [J/mol.K]
    __R = R_CONST_J__molK
    # temperature [K]
    __T_Ref = TEMPERATURE_REF_K
    # pressure [bar]
    __P_Ref = PRESSURE_REF_Pa/1e5

    def __init__(self):
        """
        Initialize the ReactionAnalyzer with a datasource and equationsource.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        """
        # self.datasource = datasource
        # self.equationsource = equationsource

    def datasource_extractor(self,
                             datasource: Dict[str, Any],
                             component_ids: List[str],
                             property_name: str):
        """
        Extract a specific property from the datasource for a given component.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        component_ids : list[str]
            The ID of the component to extract data for.
        property_name : str
            The name of the property to extract.

        Returns
        -------
        dict
            The extracted property data.
        """
        try:
            # looping through the component_id
            for component_id in component_ids:
                # check if the component exists in the datasource
                if component_id in datasource.keys():
                    # component datasource
                    component_datasource = datasource[component_id]
                    # check if the property exists in the component datasource
                    if property_name in component_datasource.keys():
                        # Extract the property data from the datasource
                        return datasource[component_id][property_name]

            # If the property is not found in any of the components, return None
            raise ValueError(
                f"Property '{property_name}' not found for component '{component_ids}'.")
        except KeyError as e:
            raise KeyError(
                f"Property '{property_name}' not found for component '{component_id}'.") from e

    def equationsource_extractor(self,
                                 equationsource: Dict[str, Any],
                                 component_ids: List[str],
                                 equation_name: str):
        """
        Extract a specific equation name from the equationsource for a given component.

        Parameters
        ----------
        equationsource : dict
            The equationsource containing the reaction equations.
        component_ids : List[str]
            The ID of the component to extract data for.
        equation_name : str
            The name of the equation to extract.


        """
        try:
            # looping through the component_id
            for component_id in component_ids:
                # check if the component exists in the equationsource
                if component_id in equationsource.keys():
                    # component equationsource
                    component_equationsource = equationsource[component_id]
                    # check if the equation name exists in the component equationsource
                    if equation_name in component_equationsource.keys():
                        # Extract the equation data from the equationsource
                        return equationsource[component_id][equation_name]

            # If the equation name is not found in any of the components, return None
            raise ValueError(
                f"Equation '{equation_name}' not found for component '{component_ids}'.")
        except KeyError as e:
            raise KeyError(
                f"Equation '{equation_name}' not found for component '{component_id}'.") from e

    def energy_analysis(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        reaction: Dict[str, Any],
        **kwargs
    ):
        '''
        Performs energy analysis of a reaction at STP

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        reaction : dict
            The reaction to be analyzed.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        data_src : dict
            energy analysis results

        Notes
        -----
        The function performs the following calculations:
        1. calculate gibbs energy of formation
        2. calculate enthalpy of formation
        3. calculate gibbs energy of reaction
        4. calculate enthalpy of reaction
        5. calculate equilibrium constant at 298.15 K and 1 bar

        - Universal gas constant is 8.314 [J/mol.K]
        - Reference temperature is 298.15 K
        - Reference pressure is 1 bar
        '''
        # SECTION: kwargs

        # NOTE: retrieve constants
        # universal gas constant [J/mol.K]
        R = self.__R
        # temperature [K]
        T = self.__T_Ref
        # pressure [bar]
        # P = self.__P_Ref

        # NOTE: thermodb components results
        # thermodb components
        data_src = {}

        # reaction results
        # reaction name
        reaction_name = reaction['name']
        # reaction
        reaction_body = reaction['reaction']

        # update
        data_src = {
            'name': reaction_name,
            'reaction': reaction_body,
            GIBBS_FREE_ENERGY_OF_FORMATION_STD: {
                'reactants': {},
                'products': {},
                'symbol': GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL,
                'unit': 'kJ/mol'
            },
            ENTHALPY_OF_FORMATION_STD: {
                'reactants': {},
                'products': {},
                'symbol': ENTHALPY_OF_FORMATION_STD_SYMBOL,
                'unit': 'kJ/mol'
            }
        }

        # SECTION: retrieve thermodynamic data
        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # molecule
            molecule_ = reactant['molecule']
            molecule_state_ = reactant['molecule_state']

            # ! gibbs energy of formation
            _dGf_IG = self.datasource_extractor(
                # type: ignore
                datasource,
                [molecule_, molecule_state_],
                'GiEnFo'
            )
            # check
            if _dGf_IG is None or _dGf_IG == 'None':
                raise ValueError(
                    f"Failed to extract Gibbs energy of formation for {molecule_}.")

            # ? add
            data_src[GIBBS_FREE_ENERGY_OF_FORMATION_STD]['reactants'][molecule_state_] = {
                'value': float(_dGf_IG['value']),
                'unit': _dGf_IG['unit']
            }

            # ! enthalpy of formation
            _dHf_IG = self.datasource_extractor(
                # type: ignore
                datasource,
                [molecule_, molecule_state_],
                'EnFo'
            )

            # check
            if _dHf_IG is None or _dHf_IG == 'None':
                raise ValueError(
                    f"Failed to extract Enthalpy of formation for {molecule_}.")

            # ? add
            data_src[ENTHALPY_OF_FORMATION_STD]['reactants'][molecule_state_] = {
                'value': float(_dHf_IG['value']),
                'unit': _dHf_IG['unit']
            }

        # NOTE: looping through products
        for product in reaction['products']:
            # molecule
            molecule_ = product['molecule']
            molecule_state_ = product['molecule_state']

            # ! gibbs energy of formation
            _dGf_IG = self.datasource_extractor(
                # type: ignore
                datasource,
                [molecule_, molecule_state_],
                'GiEnFo'
            )

            # check
            if _dGf_IG is None or _dGf_IG == 'None':
                raise ValueError(
                    f"Failed to extract Gibbs energy of formation for {molecule_}.")

            # ? add
            data_src[GIBBS_FREE_ENERGY_OF_FORMATION_STD]['products'][molecule_state_] = {
                'value': float(_dGf_IG['value']),
                'unit': _dGf_IG['unit']
            }

            # ! enthalpy of formation
            _dHf_IG = self.datasource_extractor(
                # type: ignore
                datasource,
                [molecule_, molecule_state_],
                'EnFo'
            )

            # check
            if _dHf_IG is None or _dHf_IG == 'None':
                raise ValueError(
                    f"Failed to extract Enthalpy of formation for {molecule_}.")

            # ? add
            data_src[ENTHALPY_OF_FORMATION_STD]['products'][molecule_state_] = {
                'value': float(_dHf_IG['value']),
                'unit': _dHf_IG['unit']
            }

        # SECTION: calculate gibbs energy of reaction
        # # NOTE: gibbs energy of reaction
        gibbs_energy_of_reaction_item = 0
        enthalpy_of_reaction_item = 0

        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # molecule
            molecule_ = reactant['molecule']
            # molecule state
            molecule_state_ = reactant['molecule_state']

            # ! calculate gibbs energy of reaction
            val_0 = data_src[GIBBS_FREE_ENERGY_OF_FORMATION_STD]['reactants'][molecule_state_]['value']
            gibbs_energy_of_reaction_item -= val_0 * \
                reactant['coefficient']

            # ! calculate enthalpy of reaction
            val_1 = data_src[ENTHALPY_OF_FORMATION_STD]['reactants'][molecule_state_]['value']
            enthalpy_of_reaction_item -= val_1 * \
                reactant['coefficient']

        # NOTE: looping through products
        for product in reaction['products']:
            # molecule
            molecule_ = product['molecule']
            # molecule state
            molecule_state_ = product['molecule_state']

            # ! calculate gibbs energy of reaction
            val_0 = data_src[GIBBS_FREE_ENERGY_OF_FORMATION_STD]['products'][molecule_state_]['value']
            gibbs_energy_of_reaction_item += val_0 * \
                product['coefficient']

            # ! calculate enthalpy of reaction
            val_1 = data_src[ENTHALPY_OF_FORMATION_STD]['products'][molecule_state_]['value']
            enthalpy_of_reaction_item += val_1 * \
                product['coefficient']

        # NOTE: save results
        # Gibbs energy of reaction at 298.15 K [kJ/mol]
        data_src[GIBBS_FREE_ENERGY_OF_REACTION_STD] = {
            'value': float(gibbs_energy_of_reaction_item),
            'symbol': GIBBS_FREE_ENERGY_OF_REACTION_STD_SYMBOL,
            'unit': 'kJ/mol'
        }
        # enthalpy of reaction at 298.15 K [kJ/mol]
        data_src[ENTHALPY_OF_REACTION_STD] = {
            'value': float(enthalpy_of_reaction_item),
            'symbol': ENTHALPY_OF_REACTION_STD_SYMBOL,
            'unit': 'kJ/mol'
        }

        # NOTE: equilibrium constant at 298.15 K and 1 bar
        _val_Ka = exp(-1*gibbs_energy_of_reaction_item*1000/(R*T))
        # ? add
        data_src[EQUILIBRIUM_CONSTANT_STD] = {
            'value': float(_val_Ka),
            'symbol': EQUILIBRIUM_CONSTANT_STD_SYMBOL,
            'unit': 'dimensionless'
        }

        # res
        return data_src

    def component_energy_at_temperature(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        component_names: List[str],
        temperature: float,
        res_format: Literal[
            'symbolic', 'names'
        ] = 'names',
        **kwargs
    ):
        """
        Calculate Gibbs and enthalpy energies at a given temperature.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        component_names : List[str]
            The names of the components for which to calculate Gibbs energy and enthalpy.
        temperature : float
            The temperature at which to calculate Gibbs energy.
        res_format : str, optional
            The format of the result, by default 'names'.
        kwargs : dict
            Additional keyword arguments.
            - decimal_accuracy : int
                Set decimal accuracy.

        Returns
        -------
        dict
            A dictionary containing the Gibbs energy and enthalpy at the given
            temperature as:
            - EnFo: Enthalpy of formation at 298.15 K.
            - GiEnFo: Gibbs energy of formation at 298.15 K.
            - GiEn_T: Gibbs energy at temperature T.
            - En_T: Enthalpy at temperature T.

        Notes
        -----
        The function calculates the Gibbs energy at a given temperature using the following equation:

        #### (ΔG° / RT) = (ΔG°₀ - ΔH°₀) / (R * T₀) + ΔH°₀ / (R * T) + (1 / T) * ∫(T₀ to T) [ (ΔC°p / R) dT ] - ∫(T₀ to T) [ (ΔC°p / R) * (dT / T) ]

        Where:
        - ΔG° : Standard Gibbs energy change at temperature T
        - ΔG°₀ : Standard Gibbs energy change at reference temperature T₀
        - ΔH°₀ : Standard enthalpy change at reference temperature T₀
        - ΔC°p : Standard heat capacity change (Cp_products - Cp_reactants)
        - R : Gas constant
        - T : Temperature of interest
        - T₀ : Reference temperature

        The results are returned in a dictionary with the following keys:

        EnFo : dict
            Enthalpy of formation at 298.15 K.
            - value : float
                Enthalpy of formation value [J/mol].
            - unit : str
                Unit of enthalpy of formation.

        GiEnFo : dict
            Gibbs energy of formation at 298.15 K.
            - value : float
                Gibbs energy of formation value [J/mol].
            - unit : str
                Unit of Gibbs energy of formation.

        GiEn_T : dict
            Gibbs energy at temperature T.
            - value : float
                Gibbs energy value [J/mol].
            - unit : str
                Unit of Gibbs energy.

        En_T : dict
            Enthalpy at temperature T.
            - value : float
                Enthalpy value [J/mol].
            - unit : str
                Unit of enthalpy.

        Reference
        ----------
        - Introduction to Chemical Engineering Thermodynamics, 9th Edition, (page 544)
        """
        try:
            # NOTE: retrieve constants
            # universal gas constant [J/mol.K]
            R = self.__R
            # temperature [K]
            T_ref = self.__T_Ref
            # pressure [bar]
            # P_ref = self.__P_Ref

            # SECTION: set
            # temperature [K]
            T = temperature

            # kwargs
            # NOTE: decimal accuracy
            # decimal_accuracy = kwargs.get('decimal_accuracy', 5)

            # SECTION: Gibbs energy at temperature
            # NOTE: enthalpy of formation at 298.15 K [kJ/mol]
            EnFo_src = self.datasource_extractor(
                # type: ignore
                datasource,
                component_names,
                'EnFo'
            )

            # check
            if EnFo_src is None or EnFo_src == 'None':
                raise ValueError(
                    f"Failed to extract Enthalpy of formation for {component_names}.")

            # set values
            EnFo_val = float(EnFo_src['value'])
            # TODO: convert to [J/mol]
            EnFo_unit = EnFo_src['unit']
            # to [J/mol]
            # EnFo = EnFo_val*1e3
            EnFo = pycuc.to(EnFo_val, f"{EnFo_unit} => J/mol")

            # NOTE: Gibbs free energy of formation at 298.15 K [kJ/mol]
            GiEnFo_src = self.datasource_extractor(
                # type: ignore
                datasource,
                component_names,
                'GiEnFo'
            )

            # check
            if GiEnFo_src is None or GiEnFo_src == 'None':
                raise ValueError(
                    f"Failed to extract Gibbs energy of formation for {component_names}.")

            # set values
            GiEnFo_val = float(GiEnFo_src['value'])
            GiEnFo_unit = GiEnFo_src['unit']
            # to [J/mol]
            # GiEnFo = GiEnFo_val*1e3
            GiEnFo = pycuc.to(GiEnFo_val, f"{GiEnFo_unit} => J/mol")

            # NOTE: set equation
            # ! extract Cp equation
            _eq = self.equationsource_extractor(
                equationsource,
                component_names,
                'Cp'
            )

            # check
            if _eq is None or _eq == 'None':
                raise ValueError(
                    f"Failed to extract Cp equation for {component_names}."
                )

            # check format
            if not isinstance(_eq, TableEquation):
                raise ValueError(
                    f"Invalid Cp equation format for {component_names}."
                )

            # NOTE: integral [Cp/RT]
            # unit: [dimensionless]
            # init unit
            unit_ = None
            # scipy integrate method

            def integrand_0(T):
                '''
                Dimensionless integrand for Cp/RT calculation.
                '''
                # modify the nonlocal variable
                nonlocal unit_
                res_ = _eq.cal(T=T)
                cal_ = res_.get('value', None)
                unit_ = res_.get('unit', None)
                if cal_ is None:
                    raise ValueError(
                        f"Failed to calculate Cp for {component_names} at T={T} K.")

                if not isinstance(cal_, str) and not isinstance(cal_, float):
                    raise ValueError(
                        f"Invalid Cp value for {component_names} at T={T} K: {cal_}")

                if unit_ is None:
                    raise ValueError(
                        f"Failed to get unit for Cp of {component_names} at T={T} K.")

                res = cal_/(T*R)
                return res

            # ! check custom
            custom_integral = _eq.custom_integral
            if custom_integral.get('Cp/RT', None) is not None:
                # calc
                # method 1
                _eq_Cp_integral_Cp__RT = _eq.cal_custom_integral(
                    'Cp/RT',
                    T1=T_ref,
                    T2=T
                )
            else:
                # calc
                # method 2
                _eq_Cp_integral_Cp__RT, _ = integrate.quad(
                    integrand_0,
                    T_ref,
                    T
                )
                # check
                if unit_ is None:
                    raise ValueError(
                        f"equation {component_names} does not have a unit for Cp integral.")

                # TODO: FINALLY convert Cp equation to [J/mol.K]
                _eq_Cp_integral_Cp__RT = pycuc.to(_eq_Cp_integral_Cp__RT,
                                                  f"{unit_} => J/mol.K")

            # check
            if not _eq_Cp_integral_Cp__RT:
                raise ValueError(
                    f"Failed to calculate Cp integral for {component_names}.")

            # NOTE: Cp integral
            # unit: [J/mol]
            # init unit_
            unit_ = None

            # scipy integrate method
            def integrand_1(T):
                '''
                Integral of Cp equation, unit: [J/mol]
                '''
                nonlocal unit_
                res_ = _eq.cal(T=T)
                cal_ = res_.get('value', None)
                unit_ = res_.get('unit', None)
                if cal_ is None:
                    raise ValueError(
                        f"Failed to calculate Cp for {component_names} at T={T} K.")

                if not isinstance(cal_, str) and not isinstance(cal_, float):
                    raise ValueError(
                        f"Invalid Cp value for {component_names} at T={T} K: {cal_}")

                if unit_ is None:
                    raise ValueError(
                        f"Failed to get unit for Cp of {component_names} at T={T} K.")

                return cal_

            # ! check integral
            function_integral = _eq.body_integral

            # check
            if function_integral:
                # calc
                # method 1
                _eq_Cp_integral = _eq.cal_integral(
                    T1=T_ref,
                    T2=T
                )
            else:
                # calc
                # method 2
                _eq_Cp_integral, _ = integrate.quad(
                    integrand_1,
                    T_ref,
                    T
                )

                # check
                if unit_ is None:
                    raise ValueError(
                        f"equation {component_names} does not have a unit for Cp integral.")

                # TODO: FINALLY convert Cp equation to [J/mol.K] (after integration means the unit is [J/mol])
                _eq_Cp_integral = pycuc.to(
                    _eq_Cp_integral, f"{unit_} => J/mol.K")

            # check
            if not _eq_Cp_integral:
                raise ValueError(
                    f"Failed to calculate Cp integral for {component_names}.")

            # ! enthalpy of formation at T [J/mol]
            En_T = float(EnFo + (_eq_Cp_integral))

            # ! Gibbs free energy of formation at T [J/mol]
            # A [J/mol]
            A = (GiEnFo - EnFo)/(R*T_ref)
            # B
            B = EnFo/(R*T)
            # C
            C = (1/T)*_eq_Cp_integral/R
            # D
            D = _eq_Cp_integral_Cp__RT
            # E
            E = A + B + C - D
            # at T [J/mol]
            GiEn_T = float(E*R*T)

            # SECTION: set results
            # check format
            if res_format == "names":
                ENTHALPY_OF_FORMATION_STD_RES = ENTHALPY_OF_FORMATION_STD
                GIBBS_FREE_ENERGY_OF_FORMATION_STD_RES = GIBBS_FREE_ENERGY_OF_FORMATION_STD
                ENTHALPY_OF_FORMATION_T_RES = ENTHALPY_OF_FORMATION_T
                GIBBS_FREE_ENERGY_OF_FORMATION_T_RES = GIBBS_FREE_ENERGY_OF_FORMATION_T
            elif res_format == "symbolic":
                ENTHALPY_OF_FORMATION_STD_RES = ENTHALPY_OF_FORMATION_STD_SYMBOL
                GIBBS_FREE_ENERGY_OF_FORMATION_STD_RES = GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL
                ENTHALPY_OF_FORMATION_T_RES = ENTHALPY_OF_FORMATION_T_SYMBOL
                GIBBS_FREE_ENERGY_OF_FORMATION_T_RES = GIBBS_FREE_ENERGY_OF_FORMATION_T_SYMBOL
            else:
                raise ValueError(
                    f"Invalid result format: {res_format}. Use 'names' or 'symbolic'.")

            return {
                ENTHALPY_OF_FORMATION_STD_RES: {
                    'value': EnFo,
                    'unit': 'J/mol',
                    'symbol': ENTHALPY_OF_FORMATION_STD_SYMBOL,
                    'name': ENTHALPY_OF_FORMATION_STD
                },
                GIBBS_FREE_ENERGY_OF_FORMATION_STD_RES: {
                    'value': GiEnFo,
                    'unit': 'J/mol',
                    'symbol': GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL,
                    'name': GIBBS_FREE_ENERGY_OF_FORMATION_STD
                },
                ENTHALPY_OF_FORMATION_T_RES: {
                    'value': En_T,
                    'unit': 'J/mol',
                    'symbol': ENTHALPY_OF_FORMATION_T_SYMBOL,
                    'name': ENTHALPY_OF_FORMATION_T
                },
                GIBBS_FREE_ENERGY_OF_FORMATION_T_RES: {
                    'value': GiEn_T,
                    'unit': 'J/mol',
                    'symbol': GIBBS_FREE_ENERGY_OF_FORMATION_T_SYMBOL,
                    'name': GIBBS_FREE_ENERGY_OF_FORMATION_T
                },
            }
        except Exception as e:
            raise Exception(
                f"Error in ReactionAnalyzer.Gibbs_energy_at_temperature(): {str(e)}") from e

    def reaction_energy_analysis(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        temperature: float,
        reaction: dict,
        **kwargs
    ):
        '''
        Calculates change in Gibbs free energy and enthalpy of a reaction at different temperatures.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        temperature : float
            The temperature [K] at which to calculate Gibbs energy.
        reaction : dict
            The reaction to be analyzed.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        list
            change in Gibbs free energy and enthalpy of a reaction at different temperatures.
        '''
        # SECTION: kwargs

        # NOTE: retrieve constants
        # universal gas constant [J/mol.K]
        R = self.__R
        # temperature [K]
        # T_ref = self.__T_Ref
        # pressure [bar]
        # P_ref = self.__P_Ref

        # set temperature [K]
        T = temperature

        # SECTION: thermodb components results
        thermodb = {}

        # reaction name
        reaction_name = reaction['name']
        # reaction body
        reaction_body = reaction['reaction']

        # init
        thermodb = {
            'reaction': reaction_name,
            'reaction-body': reaction_body,
        }

        # ? update
        thermodb['temperature'] = {
            'value': float(T),
            'symbol': 'T',
            'unit': 'K'
        }
        thermodb['parms'] = {
            'reactants': {},
            'products': {}
        }

        # SECTION: reactant energy analysis
        # looping through reactants
        for reactant in reaction['reactants']:
            # molecule name
            molecule_ = reactant['molecule']
            # molecule state
            molecule_state_ = reactant['molecule_state']

            # calculate Gibbs energy and enthalpy at T
            res__ = self.component_energy_at_temperature(
                datasource,
                equationsource,
                [molecule_, molecule_state_],
                T
            )

            # save
            thermodb['parms']['reactants'][molecule_state_] = res__

        # SECTION: product energy analysis
        # looping through products
        for product in reaction['products']:
            # molecule name
            molecule_ = product['molecule']
            # molecule state
            molecule_state_ = product['molecule_state']

            # calculate Gibbs energy and enthalpy at T
            res__ = self.component_energy_at_temperature(
                datasource,
                equationsource,
                [molecule_, molecule_state_],
                T
            )

            # save
            thermodb['parms']['products'][molecule_state_] = res__

        # SECTION: calculate Gibbs energy of reaction
        # overall energy analysis
        _val_dGrxn_T = 0
        _val_dHrxn_T = 0

        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # reactant name
            reactant_name = reactant['molecule_state']
            # src
            src_ = thermodb['parms']['reactants'][reactant_name]

            # dGrxn at T
            _val_dGrxn_T -= src_[GIBBS_FREE_ENERGY_OF_FORMATION_T]['value'] * \
                reactant['coefficient']
            # dHrxn at T
            _val_dHrxn_T -= src_[ENTHALPY_OF_FORMATION_T]['value'] * \
                reactant['coefficient']

        # NOTE: looping through products
        for product in reaction['products']:
            # product name
            product_name = product['molecule_state']
            # src
            src_ = thermodb['parms']['products'][product_name]

            # dGrxn at T
            _val_dGrxn_T += src_[GIBBS_FREE_ENERGY_OF_FORMATION_T]['value'] * \
                product['coefficient']

            # dHrxn at T
            _val_dHrxn_T += src_[ENTHALPY_OF_FORMATION_T]['value'] * \
                product['coefficient']

        # NOTE: save
        thermodb[GIBBS_FREE_ENERGY_OF_REACTION_T] = {
            'value': float(_val_dGrxn_T),
            'symbol': GIBBS_FREE_ENERGY_OF_REACTION_T_SYMBOL,
            'unit': 'J/mol'
        }
        thermodb[ENTHALPY_OF_REACTION_T] = {
            'value': float(_val_dHrxn_T),
            'symbol': ENTHALPY_OF_REACTION_T_SYMBOL,
            'unit': 'J/mol'
        }

        # res
        return thermodb

    def vh(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        temperature: float,
        reaction: dict,
        **kwargs
    ):
        '''
        Calculates change in Gibbs free energy of a reaction at different temperatures using the Van't Hoff equation.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        temperature : float
            The temperature [K] at which to calculate Gibbs energy.
        reaction : dict
            The reaction to be analyzed.
        kwargs : dict
            Additional keyword arguments.


        Returns
        -------
        list
            change in Gibbs free energy of a reaction at different temperatures
        '''
        try:
            # SECTION: kwargs

            # NOTE: retrieve constants
            # universal gas constant [J/mol.K]
            R = self.__R
            # temperature [K]
            # T_ref = self.__T_Ref
            # pressure [bar]
            # P_ref = self.__P_Ref

            # set temperature [K]
            T = temperature

            # SECTION: reactant energy analysis
            res_ = self.reaction_energy_analysis(
                datasource,
                equationsource,
                temperature,
                reaction
            )

            # NOTE: Gibbs energy of reaction at T [J/mol]
            _val_dGrxn_T = res_[GIBBS_FREE_ENERGY_OF_REACTION_T]['value']

            # SECTION: equilibrium constant
            Ka = exp(-1*_val_dGrxn_T/(R*T))
            # ? save
            return {
                'value': float(Ka),
                "symbol": EQUILIBRIUM_CONSTANT_SYMBOL,
                'unit': 'dimensionless',
                'temperature': {
                    'value': float(T),
                    'symbol': 'T',
                    'unit': 'K'
                },
                'reaction': {
                    'name': reaction['name'],
                    'reaction': reaction['reaction']
                },
            }
        except Exception as e:
            raise Exception(
                f"Error in ReactionAnalyzer.vh(): {str(e)}") from e

    def vh_shortcut(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        temperature: float,
        enthalpy_of_reaction_std: float,
        equilibrium_constant_std: float,
        reaction_name: str,
        reaction_body: str,
        **kwargs
    ):
        """
        Shortcut for Van't Hoff equation.

        Parameters
        ----------
        datasource : dict
            The datasource containing the thermodynamic data.
        equationsource : dict
            The equationsource containing the reaction equations.
        temperature : float
            The temperature [K] at which to calculate Gibbs energy.
        enthalpy_of_reaction_std : float
            Enthalpy of reaction at standard conditions [kJ/mol].
        equilibrium_constant_std : float
            Equilibrium constant at standard conditions [dimensionless].
        reaction_name : str
            Name of the reaction.
        reaction_body : str
            Reaction body.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        float
            Equilibrium constant at the given temperature [dimensionless].
        """
        try:
            # NOTE:convert energy to [J/mol]
            A = (-1*enthalpy_of_reaction_std*1000/self.__R) * \
                (1/temperature - 1/self.__T_Ref)
            Ka = equilibrium_constant_std*exp(A)

            # res
            return {
                'value': float(Ka),
                'symbol': EQUILIBRIUM_CONSTANT_SYMBOL,
                'unit': 'dimensionless',
                'temperature': {
                    'value': float(temperature),
                    'symbol': 'T',
                    'unit': 'K'
                },
                'reaction': {
                    'name': reaction_name,
                    'reaction': reaction_body
                },
            }

        except Exception as e:
            raise Exception(
                f"Error in ReactionAnalyzer.vh_shortcut(): {str(e)}") from e

    @staticmethod
    def norm_mole_fraction(mole_fraction: Dict[str, float | int]):
        """
        Normalize mole fractions (Xf) to sum to 1.

        Parameters
        ----------
        mole_fraction : dict
            Dictionary with species as keys and mole fractions as values.

        Returns
        -------
        dict
            Dictionary with species as keys and normalized mole fractions as values.
        """
        # Calculate total mole fraction
        total_mole_fraction = sum(mole_fraction.values())

        # Check if total mole fraction is zero to avoid division by zero
        if total_mole_fraction == 0:
            raise ValueError("Total mole fraction is zero, cannot normalize.")

        # Normalize mole fraction
        normalized_mole_fraction = {key: value / total_mole_fraction for key,
                                    value in mole_fraction.items()}

        return normalized_mole_fraction

    @staticmethod
    def cal_mole_fraction(initial_moles: Dict[str, float | int]):
        """
        Calculate mole fractions (Xf) from initial moles.

        Parameters
        ----------
        initial_moles : dict
            Dictionary with species as keys and initial moles as values.

        Returns
        -------
        mole_fraction : dict
            Dictionary with species as keys and mole fractions as values.
        total_mole_fraction : float
            Total mole fraction (should be 1.0).
        """
        # Calculate total mole
        total_mole = sum(initial_moles.values())

        # check if total mole is zero to avoid division by zero
        if total_mole == 0:
            raise ValueError(
                "Total moles are zero, cannot calculate mole fraction.")

        # Calculate mole fraction
        mole_fraction = {key: value / total_mole for key,
                         value in initial_moles.items()}

        return mole_fraction, total_mole

    @staticmethod
    def set_initial_mole(
        initial_mole: Dict[str, float | int],
        minimum_mole: float = 1e-5
    ):
        """
        Set initial moles (Xi) of each species in a reaction system.

        Parameters
        ----------
        initial_mole : dict
            Dictionary with species as keys and initial moles as values.

        Returns
        -------

        """
        try:
            # set minimum for zero mole
            initial_mole = {
                key: value if value >= 0 else minimum_mole for key,
                value in initial_mole.items()
            }

            return initial_mole
        except Exception as e:
            raise Exception(f"Failed to set initial moles: {str(e)}") from e

    @staticmethod
    def cal_mole(
        initial_mole_fraction: Dict[str, float | int],
        mole_basis: float = 1.0
    ):
        """
        Calculate the initial moles (Xi) of each species in a reaction system.

        Parameters
        ----------
        initial_mole_fraction : dict
            Dictionary with species as keys and initial mole fractions as values.
        mole_basis : float, optional
            The basis for calculating moles, by default 1.0.

        Returns
        -------
        dict
            Dictionary with species as keys and initial moles as values.
        """
        # Calculate total mole
        total_mole = sum(initial_mole_fraction.values()) * mole_basis

        # Calculate initial moles
        initial_moles = {key: value * total_mole for key,
                         value in initial_mole_fraction.items()}

        return initial_moles

    @staticmethod
    def set_stream(
        component_dict: Dict[str, int | float],
        mole: Dict[str, int | float],
        mole_fraction: Dict[str, int | float]
    ):
        """
        Set the stream of components in a reaction system.

        Parameters
        ----------
        component_dict : dict
            Dictionary with species as keys and their properties as values.
        mole : dict
            Dictionary with species as keys and their moles as values.
        mole_fraction : dict
            Dictionary with species as keys and their mole fractions as values.

        Returns
        -------
        tuple
            Tuple containing:
            - mole_comp_std : dict
                Dictionary with species as keys and their moles as values.
            - mole_fraction_comp_std : dict
                Dictionary with species as keys and their mole fractions as values.
            - mole_std : list
                List of moles of each species.
            - mole_fraction_std : list
                List of mole fractions of each species.
        """
        try:
            # NOTE: check if mole and mole_fraction are empty
            if not mole and not mole_fraction:
                raise ValueError("Both mole and mole_fraction are empty.")

            # NOTE: standardize mole and mole_fraction
            mole_std = []
            mole_comp_std = {}
            mole_fraction_std = []
            mole_fraction_comp_std = {}

            # NOTE: Set the stream of components
            # loop through the component_dict (id)
            for key, value in component_dict.items():
                # retrieve mole and mole_fraction
                _mole = mole[key]
                _mole_fraction = mole_fraction[key]

                # NOTE: standardize
                # ? mole
                mole_comp_std[key] = _mole
                mole_std.append(_mole)
                # ? mole fraction
                mole_fraction_comp_std[key] = _mole_fraction
                mole_fraction_std.append(_mole_fraction)

            # return
            return (mole_comp_std,
                    mole_fraction_comp_std,
                    mole_std,
                    mole_fraction_std)
        except Exception as e:
            raise Exception(f"Failed to set stream: {str(e)}") from e

    @staticmethod
    def set_phase_stream(
        initial_mole: Dict[str, float | int],
        initial_mole_fraction: Dict[str, float | int],
        phase_contents: Dict[str, List[str]]
    ):
        '''
        Set new mole and mole fraction for each phase in a reaction system.

        Parameters
        ----------
        initial_mole : dict
            Initial moles of each component.
        initial_mole_fraction : dict
            Initial mole fractions of each component.
        phase_contents : dict
            Dictionary with phase names as keys and lists of component names as values.

        '''
        try:

            # NOTE: phase stream
            phase_stream = {
                'g': {},
                'l': {},
                's': {},
                'aq': {}
            }

            # SECTION: separate components by phase
            for phase, components in phase_contents.items():

                # NOTE: check phase components exists
                if len(components) == 0:
                    continue

                # NOTE: loop through the components in the phase
                for component in components:
                    # check if component exists in initial_mole and initial_mole_fraction
                    if (component not in initial_mole or
                            component not in initial_mole_fraction):
                        raise ValueError(
                            f"Component {component} not found in initial moles or mole fractions.")

                    # set the mole and mole fraction for the phase
                    phase_stream[phase][component] = {
                        'mole': initial_mole[component],
                        'mole_fraction': initial_mole_fraction[component]
                    }

            # SECTION: calculate new moles and mole fractions for each phase
            # loop through each phase and its components
            for phase, components in phase_stream.items():

                # NOTE: check if components are empty
                if len(components) == 0:
                    continue

                # NOTE: calculate total mole for the phase
                total_mole = sum(comp['mole'] for comp in components.values())

                # check if total mole is zero to avoid division by zero
                if total_mole == 0:
                    raise ValueError(
                        f"Total moles for phase {phase} are zero, cannot calculate mole fractions.")

                # calculate new moles and mole fractions for the phase
                for component, values in components.items():
                    new_mole_ = values['mole']
                    new_mole_fraction_ = values['mole'] / total_mole

                    # update the new_mole and new_mole_fraction dictionaries
                    phase_stream[phase][component]['phase_mole'] = new_mole_
                    phase_stream[phase][component]['phase_mole_fraction'] = new_mole_fraction_

            # return the new moles and mole fractions
            return phase_stream
        except Exception as e:
            raise Exception(f"Failed to set phase stream: {str(e)}") from e

    @staticmethod
    def cal_conversion(
            initial_mole: Dict[str, float | int],
            final_mole: Dict[str, float | int],
            components: List[str]
    ):
        '''
        Calculate conversion

        Parameters
        ----------
        initial_mole : dict
            Initial moles of each component.
        final_mole : dict
            Final moles of each component.
        components : list
            List of components to calculate conversion for.

        Returns
        -------

        '''
        try:
            # check if N0s and Nfs are empty
            if not initial_mole and not final_mole:
                raise ValueError("Both N0s and Nfs are empty.")

            # check if components is empty
            if not components:
                raise ValueError("Components list is empty.")

            # NOTE: calculate conversion
            # results
            conversion = {}

            # loop through the component_item
            for component in components:
                # initial mole
                _N0i = initial_mole[component]
                # final mole
                _Nfi = final_mole[component]
                # conversion
                conversion[component] = {
                    'value': ((_N0i - _Nfi)/_N0i)*100,
                    'unit': '%'
                }

            return conversion
        except Exception as e:
            raise Exception(f"Failed to calculate conversion: {str(e)}") from e

    def calc_actual_gibbs_energy_of_reaction(
        self,
        reaction: dict,
        chemical_potential: Dict[str, float | int],
        temperature: Temperature,
        pressure: Pressure,
        **kwargs
    ):
        '''
        Calculates the actual Gibbs energy of reaction from the chemical potential.

        Parameters
        ----------
        reaction : dict
            The reaction to be analyzed.
        chemical_potential : dict
            Dictionary containing the chemical potential of each component in the reaction.
        temperature : Temperature
            The temperature at which to calculate the Gibbs energy of reaction.
        pressure : Pressure
            The pressure at which to calculate the Gibbs energy of reaction.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        list
            change in Gibbs free energy and enthalpy of a reaction at different temperatures.
        '''
        # SECTION: kwargs

        # NOTE: retrieve constants
        # universal gas constant [J/mol.K]
        R = self.__R
        # temperature [K]
        # T_ref = self.__T_Ref
        # pressure [bar]
        # P_ref = self.__P_Ref

        # set temperature [K]
        T = temperature

        # SECTION: thermodb components results
        thermodb = {}

        # reaction name
        reaction_name = reaction['name']
        # reaction body
        reaction_body = reaction['reaction']

        # init
        thermodb = {
            'reaction': reaction_name,
            'reaction-body': reaction_body,
        }

        # ? update
        thermodb['temperature'] = {
            'value': temperature['value'],
            'symbol': 'T',
            'unit': temperature['unit']
        }
        thermodb['pressure'] = {
            'value': pressure['value'],
            'symbol': 'P',
            'unit': pressure['unit']
        }
        thermodb['parms'] = {
            'reactants': {},
            'products': {}
        }

        # SECTION: reactant energy analysis
        # looping through reactants
        for reactant in reaction['reactants']:
            # molecule name
            molecule_ = reactant['molecule']
            # molecule state
            molecule_state_ = reactant['molecule_state']

            # retrieve chemical potential
            res__ = chemical_potential.get(
                molecule_state_, None)

            # check if chemical potential is None
            if res__ is None:
                raise ValueError(
                    f"Chemical potential for {molecule_state_} not found in the provided data.")

            # save
            thermodb['parms']['reactants'][molecule_state_] = {
                CHEMICAL_POTENTIAL_MIXTURE_T_P: {
                    'value': res__,
                    'symbol': CHEMICAL_POTENTIAL_MIXTURE_T_P_SYMBOL,
                    'unit': 'J/mol'
                }
            }

        # SECTION: product energy analysis
        # looping through products
        for product in reaction['products']:
            # molecule name
            molecule_ = product['molecule']
            # molecule state
            molecule_state_ = product['molecule_state']

            # calculate chemical potential
            res__ = chemical_potential.get(
                molecule_state_, None)

            # check if chemical potential is None
            if res__ is None:
                raise ValueError(
                    f"Chemical potential for {molecule_state_} not found in the provided data.")

            # save
            thermodb['parms']['products'][molecule_state_] = {
                CHEMICAL_POTENTIAL_MIXTURE_T_P: {
                    'value': res__,
                    'symbol': CHEMICAL_POTENTIAL_MIXTURE_T_P_SYMBOL,
                    'unit': 'J/mol'
                }
            }

        # SECTION: calculate actual Gibbs energy of reaction
        # overall energy analysis [J/mol]
        _val_dGrxn_T = 0

        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # reactant name
            reactant_name = reactant['molecule_state']
            # src
            src_ = thermodb['parms']['reactants'][reactant_name]

            # dGrxn at T
            _val_dGrxn_T -= src_[CHEMICAL_POTENTIAL_MIXTURE_T_P]['value'] * \
                reactant['coefficient']

        # NOTE: looping through products
        for product in reaction['products']:
            # product name
            product_name = product['molecule_state']
            # src
            src_ = thermodb['parms']['products'][product_name]

            # dGrxn at T
            _val_dGrxn_T += src_[CHEMICAL_POTENTIAL_MIXTURE_T_P]['value'] * \
                product['coefficient']

        # NOTE: save
        thermodb[ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION] = {
            'value': float(_val_dGrxn_T),
            'symbol': ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION_SYMBOL,
            'unit': 'J/mol'
        }

        # res
        return thermodb
