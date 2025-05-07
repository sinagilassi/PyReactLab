# import libs
from typing import Dict, Any, List, Literal, Optional
from math import exp
#  local
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K,
    EnFo_IG, EnFo_LIQ, GiEnFo_IG, GiEnFo_LIQ
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

    def energy_analysis(self,
                        datasource: Dict[str, Any],
                        equationsource: Dict[str, Any],
                        reaction: Dict[str, Any],
                        **kwargs):
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
            - decimal_accuracy : int

        Returns
        -------
        thermodb_component : dict
            thermodb component data

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
        # NOTE: decimal accuracy
        decimal_accuracy = kwargs.get('decimal_accuracy', 5)

        # NOTE: retrieve constants
        # universal gas constant [J/mol.K]
        R = self.__R
        # temperature [K]
        T = self.__T_Ref
        # pressure [bar]
        P = self.__P_Ref

        # NOTE: thermodb components results
        # thermodb components
        thermodb_component = {}

        # reaction results
        # reaction name
        reaction_name = reaction['name']
        # reaction
        reaction_body = reaction['reaction']

        # update
        thermodb_component = {
            'name': reaction_name,
            'reaction': reaction_body,
            'dGf_IG': {
                'reactants': {},
                'products': {}
            },
            'dHf_IG': {
                'reactants': {},
                'products': {}
            },
            'dGrxn_298': 0,
            'dHrxn_298': 0,
            'Ka': 0
        }

        # SECTION: retrieve thermodynamic data
        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # molecule
            molecule_ = reactant['molecule']
            # ! gibbs energy of formation
            _dGf_IG = datasource[molecule_][GiEnFo_IG]
            _val_dGf_IG = _dGf_IG['value']
            # _val_dGf_IG = thermodb[reactant['molecule']].check_property(
            #     'GENERAL').get_property('dGf_IG')['value']
            # ! enthalpy of formation
            _dHf_IG = datasource[molecule_][EnFo_IG]
            _val_dHf_IG = _dHf_IG['value']
            # _val_dHf_IG = thermodb[reactant['molecule']].check_property(
            #     'GENERAL').get_property('dHf_IG')['value']

            # save
            thermodb_component['dGf_IG']['reactants'][molecule_] = float(
                _val_dGf_IG)
            thermodb_component['dHf_IG']['reactants'][molecule_] = float(
                _val_dHf_IG)

        # NOTE: looping through products
        for product in reaction['products']:
            # molecule
            molecule_ = product['molecule']
            # ! gibbs energy of formation
            _dGf_IG = datasource[molecule_][GiEnFo_IG]
            _val_dGf_IG = _dGf_IG['value']
            # _val_dGf_IG = thermodb[product['molecule']].check_property(
            #     'GENERAL').get_property('dGf_IG')['value']
            # ! enthalpy of formation
            _dHf_IG = datasource[molecule_][EnFo_IG]
            _val_dHf_IG = _dHf_IG['value']
            # _val_dHf_IG = thermodb[product['molecule']].check_property(
            #     'GENERAL').get_property('dHf_IG')['value']
            # save
            thermodb_component['dGf_IG']['products'][molecule_] = float(
                _val_dGf_IG)
            thermodb_component['dHf_IG']['products'][molecule_] = float(
                _val_dHf_IG)

        # SECTION: calculate gibbs energy of reaction
        # # NOTE: gibbs energy of reaction
        gibbs_energy_of_reaction_item = 0
        enthalpy_of_reaction_item = 0

        # NOTE: looping through reactants
        for reactant in reaction['reactants']:
            # molecule
            molecule_ = reactant['molecule']
            # ! calculate gibbs energy of reaction
            val_0 = thermodb_component['dGf_IG']['reactants'][molecule_]
            gibbs_energy_of_reaction_item -= val_0 * reactant['coefficient']

            # ! calculate enthalpy of reaction
            val_1 = thermodb_component['dHf_IG']['reactants'][molecule_]
            enthalpy_of_reaction_item -= val_1 * reactant['coefficient']

        # NOTE: looping through products
        for product in reaction['products']:
            # molecule
            molecule_ = product['molecule']
            # ! calculate gibbs energy of reaction
            val_0 = thermodb_component['dGf_IG']['products'][molecule_]
            gibbs_energy_of_reaction_item += val_0 * product['coefficient']

            # ! calculate enthalpy of reaction
            val_1 = thermodb_component['dHf_IG']['products'][molecule_]
            enthalpy_of_reaction_item += val_1 * product['coefficient']

        # NOTE: save results
        # Gibbs energy of reaction at 298.15 K [kJ/mol]
        thermodb_component['dGrxn_298'] = round(
            gibbs_energy_of_reaction_item, decimal_accuracy)
        # enthalpy of reaction at 298.15 K [kJ/mol]
        thermodb_component['dHrxn_298'] = round(
            enthalpy_of_reaction_item, decimal_accuracy)

        # NOTE: equilibrium constant at 298.15 K and 1 bar
        _val_Ka = exp(-1*gibbs_energy_of_reaction_item*1000/(R*T))
        thermodb_component['Ka'] = round(_val_Ka, decimal_accuracy)

        # res
        return thermodb_component

    def vant_hoff(self, temperatures, reaction, energy_analysis_res, thermodb, decimal_accuracy=5):
        '''
        Calculates change in Gibbs free energy of a reaction at different temperatures.

        Parameters
        ----------
        temperatures : list
            temperature range
        energy_analysis_res : dict
            energy analysis results
        reaction : dict
            reaction results
        thermodb : dict
            thermo database

        decimal_accuracy : int
            set decimal accuracy

        Returns
        -------
        list
            change in Gibbs free energy of a reaction at different temperatures
        '''
        # universal gas constant [J/mol.K]
        R = 8.314
        # Tref [K]
        Tref = 298.15
        # res
        thermodb_component = {}

        # reaction name
        reaction_name = reaction['name']
        # reaction body
        reaction_body = reaction['reaction']

        # init
        thermodb_component = {
            'reaction': reaction_name,
            'reaction-body': reaction_body,
            'Ts': {}
        }

        # looping through temperature
        for T in temperatures:
            # res
            res = {
                'T': None,
                'parms': {
                    'reactants': {},
                    'products': {}
                },
                'dGrxn_Ts': None,
                'dHrxn_Ts': None,
                'dGrxn_298': 0,
                'dHrxn_298': 0,
                'Kas': None
            }

            # temperature
            thermodb_component['Ts'][str(T)] = {}
            thermodb_component['Ts'][str(T)]['T'] = T
            thermodb_component['Ts'][str(T)]['parms'] = {
                'reactants': {},
                'products': {}
            }
            thermodb_component['Ts'][str(T)]['dGrxn_Ts'] = []
            thermodb_component['Ts'][str(T)]['dHrxn_Ts'] = []
            thermodb_component['Ts'][str(T)]['dGrxn_298'] = 0
            thermodb_component['Ts'][str(T)]['dHrxn_298'] = 0
            thermodb_component['Ts'][str(T)]['Kas'] = []

            # looping through reactants
            for reactant in reaction['reactants']:

                # symbol
                reactant_symbol = reactant['molecule']

                # enthalpy of formation at 298.15 K [J/mol]
                dHf_IG = float(
                    energy_analysis_res['dHf_IG']['reactants'][reactant_symbol])
                dHf_IG = dHf_IG*1e3
                # Gibbs free energy of formation at 298.15 K [J/mol]
                dGf_IG = float(
                    energy_analysis_res['dGf_IG']['reactants'][reactant_symbol])
                dGf_IG = dGf_IG*1e3

                # integral Cp/RT
                _eq_Cp_integral_Cp__RT = thermodb[reactant['molecule']].check_function(
                    'HEAT-CAPACITY').cal_custom_integral('Cp/RT', T1=Tref, T2=T)
                # integral Cp/R
                _eq_Cp_integral_Cp__R = thermodb[reactant['molecule']].check_function(
                    'HEAT-CAPACITY').cal_custom_integral('Cp/R', T1=Tref, T2=T)
                # Cp integral
                _eq_Cp_integral = thermodb[reactant['molecule']].check_function(
                    'HEAT-CAPACITY').cal_integral(T1=Tref, T2=T)

                # enthalpy of formation at T [J/mol]
                dHf_T = dHf_IG + (_eq_Cp_integral)
                # round
                dHf_T = round(dHf_T, decimal_accuracy)

                # Gibbs free energy of formation at T [J/mol]
                # A [J/mol]
                A = (dGf_IG - dHf_IG)/(R*Tref)
                # B
                B = dHf_IG/(R*T)
                # C
                C = (1/T)*_eq_Cp_integral/R
                # D
                D = _eq_Cp_integral_Cp__RT
                # E
                E = A + B + C - D
                # at T [kJ/mol]
                dGf_T = float(E*R*T)
                # round
                dGf_T = round(dGf_T, decimal_accuracy)

                # parms
                parms = {
                    'Cp/RT': _eq_Cp_integral_Cp__RT,
                    'Cp/R': _eq_Cp_integral_Cp__R,
                    'Cp_T': _eq_Cp_integral,
                    'dHf_T': dHf_T,
                    'dGf_T': dGf_T,
                    'A': A,
                    'B': B,
                    'C': C,
                    'D': D,
                    'E': E
                }
                # save
                thermodb_component['Ts'][str(
                    T)]['parms']['reactants'][reactant['molecule']] = parms

            # looping through products
            for product in reaction['products']:

                # product symbol
                product_symbol = product['molecule']

                # enthalpy of formation at 298.15 K [J/mol]
                dHf_IG = float(
                    energy_analysis_res['dHf_IG']['products'][product_symbol])
                dHf_IG = dHf_IG*1e3
                # Gibbs free energy of formation at 298.15 K [J/mol]
                dGf_IG = float(
                    energy_analysis_res['dGf_IG']['products'][product_symbol])
                dGf_IG = dGf_IG*1e3

                # integral Cp/RT
                _eq_Cp_integral_Cp__RT = thermodb[product['molecule']].check_function(
                    'HEAT-CAPACITY').cal_custom_integral('Cp/RT', T1=Tref, T2=T)
                # integral Cp/R
                _eq_Cp_integral_Cp__R = thermodb[product['molecule']].check_function(
                    'HEAT-CAPACITY').cal_custom_integral('Cp/R', T1=Tref, T2=T)
                # Cp integral
                _eq_Cp_integral = thermodb[product['molecule']].check_function(
                    'HEAT-CAPACITY').cal_integral(T1=Tref, T2=T)

                # enthalpy of formation at T [J/mol]
                dHf_T = dHf_IG + (_eq_Cp_integral)
                # round
                dHf_T = round(dHf_T, decimal_accuracy)

                # Gibbs free energy of formation at T [J/mol]
                # A [J/mol]
                A = (dGf_IG - dHf_IG)/(R*Tref)
                # B
                B = dHf_IG/(R*T)
                # C
                C = (1/T)*_eq_Cp_integral/R
                # D
                D = _eq_Cp_integral_Cp__RT
                # E
                E = A + B + C - D
                # at T [J/mol]
                dGf_T = float(E*R*T)
                # round
                dGf_T = round(dGf_T, decimal_accuracy)

                # parms
                parms = {
                    'Cp/RT': _eq_Cp_integral_Cp__RT,
                    'Cp/R': _eq_Cp_integral_Cp__R,
                    'Cp_T': _eq_Cp_integral,
                    'dHf_T': dHf_T,
                    'dGf_T': dGf_T,
                    'A': A,
                    'B': B,
                    'C': C,
                    'D': D,
                    'E': E
                }
                # save
                thermodb_component['Ts'][str(
                    T)]['parms']['products'][product['molecule']] = parms

            # overall
            _val_dGrxn_T = 0
            _val_dHrxn_T = 0

            # looping through reactants
            for reactant in reaction['reactants']:
                # dGrxn at T
                _val_dGrxn_T -= thermodb_component['Ts'][str(
                    T)]['parms']['reactants'][reactant['molecule']]['dGf_T']*reactant['coefficient']
                # dHrxn at T
                _val_dHrxn_T -= thermodb_component['Ts'][str(
                    T)]['parms']['reactants'][reactant['molecule']]['dHf_T']*reactant['coefficient']
            # looping through products
            for product in reaction['products']:
                # dGrxn at T
                _val_dGrxn_T += thermodb_component['Ts'][str(
                    T)]['parms']['products'][product['molecule']]['dGf_T']*product['coefficient']

                # dHrxn at T
                _val_dHrxn_T += thermodb_component['Ts'][str(
                    T)]['parms']['products'][product['molecule']]['dHf_T']*product['coefficient']
            # save
            _val_dGrxn_T = round(_val_dGrxn_T, decimal_accuracy)
            _val_dHrxn_T = round(_val_dHrxn_T, decimal_accuracy)
            thermodb_component['Ts'][str(T)]['dGrxn_Ts'] = _val_dGrxn_T
            thermodb_component['Ts'][str(T)]['dHrxn_Ts'] = _val_dHrxn_T

            # equilibrium constant
            Ka = exp(-1*float(_val_dGrxn_T)/(R*T))
            # save
            thermodb_component['Ts'][str(T)]['Kas'] = float(Ka)

        # reset
        _val_dGrxn_T = 0
        _val_dHrxn_T = 0

        # res
        return thermodb_component

    def vant_hoff_T(temperature, reaction, energy_analysis_res, thermodb, decimal_accuracy=5, R=8.314, Tref=298.15):
        '''
        Calculates change in Gibbs free energy of a reaction at different temperatures.

        Parameters
        ----------
        temperature : float
            temperature set
        energy_analysis_res : dict
            energy analysis results
        reaction : dict
            reaction results
        thermodb : dict
            thermo database
        decimal_accuracy : int
            set decimal accuracy

        Returns
        -------
        list
            change in Gibbs free energy of a reaction at different temperatures
        '''
        # universal gas constant [J/mol.K]
        # R = 8.314
        # Tref [K]
        # Tref = 298.15
        # res
        thermodb_component = {}

        # reaction name
        reaction_name = reaction['name']
        # reaction body
        reaction_body = reaction['reaction']

        # init
        thermodb_component = {
            'reaction': reaction_name,
            'reaction-body': reaction_body,
            'T': {}
        }

        # set
        T = temperature

        # res
        res = {
            'T': T,
            'parms': {
                'reactants': {},
                'products': {}
            },
            'dGrxn_Ts': None,
            'dHrxn_Ts': None,
            'dGrxn_298': 0,
            'dHrxn_298': 0,
            'Kas': None
        }

        # temperature
        thermodb_component['T'][str(T)] = {}
        thermodb_component['T'][str(T)]['T'] = T
        thermodb_component['T'][str(T)]['parms'] = {
            'reactants': {},
            'products': {}
        }
        thermodb_component['T'][str(T)]['dGrxn_Ts'] = []
        thermodb_component['T'][str(T)]['dHrxn_Ts'] = []
        thermodb_component['T'][str(T)]['dGrxn_298'] = 0
        thermodb_component['T'][str(T)]['dHrxn_298'] = 0
        thermodb_component['T'][str(T)]['Kas'] = []

        # looping through reactants
        for reactant in reaction['reactants']:

            # symbol
            reactant_symbol = reactant['molecule']

            # enthalpy of formation at 298.15 K [J/mol]
            dHf_IG = float(
                energy_analysis_res['dHf_IG']['reactants'][reactant_symbol])
            dHf_IG = dHf_IG*1e3
            # Gibbs free energy of formation at 298.15 K [J/mol]
            dGf_IG = float(
                energy_analysis_res['dGf_IG']['reactants'][reactant_symbol])
            dGf_IG = dGf_IG*1e3

            # integral Cp/RT
            _eq_Cp_integral_Cp__RT = thermodb[reactant['molecule']].check_function(
                'HEAT-CAPACITY').cal_custom_integral('Cp/RT', T1=Tref, T2=T)
            # integral Cp/R
            _eq_Cp_integral_Cp__R = thermodb[reactant['molecule']].check_function(
                'HEAT-CAPACITY').cal_custom_integral('Cp/R', T1=Tref, T2=T)
            # Cp integral
            _eq_Cp_integral = thermodb[reactant['molecule']].check_function(
                'HEAT-CAPACITY').cal_integral(T1=Tref, T2=T)

            # enthalpy of formation at T [J/mol]
            dHf_T = dHf_IG + (_eq_Cp_integral)
            # round
            dHf_T = round(dHf_T, decimal_accuracy)

            # Gibbs free energy of formation at T [J/mol]
            # A [J/mol]
            A = (dGf_IG - dHf_IG)/(R*Tref)
            # B
            B = dHf_IG/(R*T)
            # C
            C = (1/T)*_eq_Cp_integral/R
            # D
            D = _eq_Cp_integral_Cp__RT
            # E
            E = A + B + C - D
            # at T [kJ/mol]
            dGf_T = float(E*R*T)
            # round
            dGf_T = round(dGf_T, decimal_accuracy)

            # parms
            parms = {
                'Cp/RT': _eq_Cp_integral_Cp__RT,
                'Cp/R': _eq_Cp_integral_Cp__R,
                'Cp_T': _eq_Cp_integral,
                'dHf_T': dHf_T,
                'dGf_T': dGf_T,
                'A': A,
                'B': B,
                'C': C,
                'D': D,
                'E': E
            }
            # save
            thermodb_component['T'][str(
                T)]['parms']['reactants'][reactant['molecule']] = parms

        # looping through products
        for product in reaction['products']:

            # product symbol
            product_symbol = product['molecule']

            # enthalpy of formation at 298.15 K [J/mol]
            dHf_IG = float(
                energy_analysis_res['dHf_IG']['products'][product_symbol])
            dHf_IG = dHf_IG*1e3
            # Gibbs free energy of formation at 298.15 K [J/mol]
            dGf_IG = float(
                energy_analysis_res['dGf_IG']['products'][product_symbol])
            dGf_IG = dGf_IG*1e3

            # integral Cp/RT
            _eq_Cp_integral_Cp__RT = thermodb[product['molecule']].check_function(
                'HEAT-CAPACITY').cal_custom_integral('Cp/RT', T1=Tref, T2=T)
            # integral Cp/R
            _eq_Cp_integral_Cp__R = thermodb[product['molecule']].check_function(
                'HEAT-CAPACITY').cal_custom_integral('Cp/R', T1=Tref, T2=T)
            # Cp integral
            _eq_Cp_integral = thermodb[product['molecule']].check_function(
                'HEAT-CAPACITY').cal_integral(T1=Tref, T2=T)

            # enthalpy of formation at T [J/mol]
            dHf_T = dHf_IG + (_eq_Cp_integral)
            # round
            dHf_T = round(dHf_T, decimal_accuracy)

            # Gibbs free energy of formation at T [J/mol]
            # A [J/mol]
            A = (dGf_IG - dHf_IG)/(R*Tref)
            # B
            B = dHf_IG/(R*T)
            # C
            C = (1/T)*_eq_Cp_integral/R
            # D
            D = _eq_Cp_integral_Cp__RT
            # E
            E = A + B + C - D
            # at T [J/mol]
            dGf_T = float(E*R*T)
            # round
            dGf_T = round(dGf_T, decimal_accuracy)

            # parms
            parms = {
                'Cp/RT': _eq_Cp_integral_Cp__RT,
                'Cp/R': _eq_Cp_integral_Cp__R,
                'Cp_T': _eq_Cp_integral,
                'dHf_T': dHf_T,
                'dGf_T': dGf_T,
                'A': A,
                'B': B,
                'C': C,
                'D': D,
                'E': E
            }
            # save
            thermodb_component['T'][str(
                T)]['parms']['products'][product['molecule']] = parms

        # overall
        _val_dGrxn_T = 0
        _val_dHrxn_T = 0

        # looping through reactants
        for reactant in reaction['reactants']:
            # dGrxn at T
            _val_dGrxn_T -= thermodb_component['T'][str(
                T)]['parms']['reactants'][reactant['molecule']]['dGf_T']*reactant['coefficient']
            # dHrxn at T
            _val_dHrxn_T -= thermodb_component['T'][str(
                T)]['parms']['reactants'][reactant['molecule']]['dHf_T']*reactant['coefficient']
        # looping through products
        for product in reaction['products']:
            # dGrxn at T
            _val_dGrxn_T += thermodb_component['T'][str(
                T)]['parms']['products'][product['molecule']]['dGf_T']*product['coefficient']

            # dHrxn at T
            _val_dHrxn_T += thermodb_component['T'][str(
                T)]['parms']['products'][product['molecule']]['dHf_T']*product['coefficient']
        # save
        _val_dGrxn_T = round(_val_dGrxn_T, decimal_accuracy)
        _val_dHrxn_T = round(_val_dHrxn_T, decimal_accuracy)
        thermodb_component['T'][str(T)]['dGrxn_Ts'] = _val_dGrxn_T
        thermodb_component['T'][str(T)]['dHrxn_Ts'] = _val_dHrxn_T

        # equilibrium constant
        Ka = exp(-1*float(_val_dGrxn_T)/(R*T))
        # save
        thermodb_component['T'][str(T)]['Kas'] = float(Ka)

        # reset
        _val_dGrxn_T = 0
        _val_dHrxn_T = 0

        # res
        return thermodb_component

    def calculate_mole_fraction(self, initial_moles):
        """
        Calculate mole fractions from initial moles.

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

        # Calculate mole fraction
        mole_fraction = {key: value / total_mole for key,
                         value in initial_moles.items()}

        # Calculate total mole fraction (verification)
        total_mole_fraction = sum(mole_fraction.values())

        return mole_fraction, total_mole_fraction
