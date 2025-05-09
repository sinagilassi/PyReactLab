# import packages/modules
from typing import Dict, List, Literal, Optional, Any
import numpy as np
from math import sqrt, pow, exp
from scipy import optimize
import pyThermoModels as ptm
# local
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K
)


class ReactionOptimizer:
    '''Reaction Optimizer class'''

    # SECTION: variables
    _system_inputs = None
    # universal gas constant [J/mol.K]
    R = R_CONST_J__molK
    # temperature [K]
    T_Ref = TEMPERATURE_REF_K
    # pressure [Pa]
    P_Ref_Pa = PRESSURE_REF_Pa

    def __init__(self,
                 datasource: Dict[str, Any],
                 equationsource: Dict[str, Any],
                 component_dict: Dict[str, float | int],
                 comp_list: List[Dict[str, float | int]],
                 reaction_analysis: Dict,
                 overall_reaction_analysis: Dict,
                 **kwargs):
        '''
        Initialize ReactionOptimizer class

        Parameters
        ----------
        datasource : dict
            datasource dictionary
        equationsource : dict
            equationsource dictionary
        component_dict : dict
            component dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        comp_list : list
            component list [{key: value}, {key: value}, ...], such as [{CO2: -1.0, H2: -3.0, CH3OH: 1.0, H2O: 1.0}, ...]
        reaction_analysis : dict
            reaction analysis result
        overall_reaction_analysis : dict
            overall reaction analysis result
        kwargs : dict
            additional parameters
        '''
        # datasource
        self.datasource = datasource
        # equationsource
        self.equationsource = equationsource
        # set component dictionary
        self.component_dict = component_dict
        # set component list
        self.comp_list = comp_list
        # set reaction analysis
        self.reaction_analysis = reaction_analysis
        # set overall reaction analysis
        self.overall_reaction_analysis = overall_reaction_analysis

        # NOTE: kwargs
        self.threshold = kwargs.get('threshold', 1e-8)

        # NOTE: init PyThermoModels
        # init
        self.eos = ptm.eos()

    def unpack_X(self, mol_data_pack: Dict[str, float | int]):
        '''
        convert X to dict {key: value} and list [value]

        Parameters
        ----------
        mol_data_pack : dict
            component mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}

        Returns
        -------
        N0s_list : list
            The list of N0s
        N0s_vector : numpy.ndarray
            The vector of N0s
        Nf : float
            The total number of moles
        '''
        try:
            # size
            size = len(mol_data_pack)

            # lists
            N0s_list = []
            N0s_vector = np.zeros(size)

            # looping through
            for i, (com_key, com_value) in enumerate(self.component_dict.items()):
                N0s_list.append(mol_data_pack[com_key])
                N0s_vector[i] = mol_data_pack[com_key]

            # final
            Nf = np.sum(N0s_vector)

            return N0s_list, N0s_vector, Nf
        except Exception as e:
            raise Exception(
                f"Error in ReactionOptimizer.unpack_X(): {str(e)}") from e

    def build_EoR(self, EoR: List[float]):
        '''
        build EoR

        Parameters
        ----------
        EoR : list
            extent of reaction list

        Returns
        -------
        comp_value_list : list
            comp_value_list
        comp_value_matrix : numpy.ndarray
            comp_value_matrix
        '''
        try:
            # comp value list
            comp_value_list = []
            # new
            comp_list_updated = []
            # copy data
            for i, item in enumerate(self.comp_list):
                _item = {}
                for key, value in item.items():
                    _item[key] = value * EoR[i]
                comp_list_updated.append(_item)

            # define a value list
            for i, item in enumerate(comp_list_updated):
                _value_list = []
                for key, value in item.items():
                    _value_list.append(value)
                comp_value_list.append(_value_list)

            # value matrix
            comp_value_matrix = np.array(comp_value_list)

            return comp_value_list, comp_value_matrix
        except Exception as e:
            raise Exception(
                f"Error in ReactionOptimizer.build_EoR(): {str(e)}") from e

    def build_final_X(self,
                      N0s_vector: np.ndarray,
                      EoR_vector: np.ndarray):
        '''
        build final X

        Parameters
        ----------
        component_dict : dict
            component_dict
        N0s_vector : numpy.ndarray
            N0s_vector

        Returns
        -------
        Xfs : dict
            Xfs
        Xfs_vector : numpy.ndarray
            Xfs_vector
        Nf : float
            Nf
        '''
        try:
            # final mole of components [mole]
            Nfs_vector = N0s_vector + EoR_vector

            # final mole fraction
            Xfs_vector = Nfs_vector/np.sum(Nfs_vector)

            # final mole
            Nf = np.sum(Nfs_vector)

            # convert to Xfs
            Xfs = {}
            for i, key in enumerate(self.component_dict.keys()):
                Xfs[key] = Xfs_vector[i]

            return Xfs, Xfs_vector, Nf, Nfs_vector
        except Exception as e:
            raise Exception(
                f"Error in ReactionOptimizer.build_final_X(): {str(e)}") from e

    def equilibrium_reaction_objective_function(self,
                                                x,
                                                N0s: Dict[str, float],
                                                P: float,
                                                T: float,
                                                Keq_dict: Dict[str, float],
                                                ):
        '''
        Equilibrium Reaction Analysis

        Parameters
        ----------
        x : list
            loop values
        N0s : dict
            initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
            parameters
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        Keq_dict : dict
            reaction equilibrium constant dictionary calculated at T as: {key: value}, such as {R1: 0.1, R2: 0.2, ...}

        Returns
        -------
        obj : float
            objective function
        '''
        try:
            # NOTE: default values
            # pressure [Pa]
            P_ref = self.P_Ref_Pa

            # NOTE: extent of reaction
            EoR = x

            # extent of reaction
            EoR = x
            # define threshold
            for i, item in enumerate(EoR):
                EoR[i] = max(self.threshold, item)

            # SECTION: mole balance
            # NOTE: initial mole
            # unpack N0s
            N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

            # NOTE: update comp_list with EoR
            # build EoR => item[key] = value * EoR[i]
            _, comp_value_matrix = self.build_EoR(EoR)

            # NOTE: extent sun [mol]
            EoR_vector = np.sum(comp_value_matrix, axis=0)

            # NOTE: build final X
            Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
                N0s_vector, EoR_vector)

            # SECTION: component state analysis

            # SECTION: fugacity coefficient
            # calculate fugacity coefficient
            # model input
            # # component list
            # comp_list = list(self.component_dict.keys())

            # # feed spec
            # feed_spec = {}
            # for comp in comp_list:
            #     feed_spec[comp] = Xfs_vector[comp_list.index(comp)]

            # # model input
            # model_input = {
            #     "eos-model": eos_model,
            #     "phase": phase,
            #     "feed-spec": feed_spec,
            #     "operating-conditions": {
            #         "pressure": [P, 'Pa'],
            #         "temperature": [T, 'K'],
            #     },
            # }

            # # calculate fugacity coefficient
            # phis = []
            # # check
            # if mode == "non-ideal":
            #     # fugacity coefficient [-]
            #     fugacity_res = fugacity_obj.cal_fugacity_coefficient(
            #         model_input)
            #     # phii
            #     _, _, _, phi_pack = fugacity_res

            #     # phis
            #     phis = {}
            #     for i, key in enumerate(component_dict.keys()):
            #         phis[key] = phi_pack['VAPOR'][key]['phi']

            # elif mode == 'ideal':
            #     phis = {}
            #     for i, key in enumerate(component_dict.keys()):
            #         phis[key] = 1
            # else:
            #     raise ValueError(
            #         "Invalid mode. Must be 'ideal' or 'non-ideal'.")

            # SECTION: equilibrium equation
            # reaction term
            reaction_terms = {}

            # NOTE: loop over reactions
            for i, reaction in enumerate(self.reaction_analysis):
                # denominator
                denominator = 1
                # numerator
                numerator = 1

                # SECTION: loop over reactants
                for item in self.reaction_analysis[reaction]['reactants']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_]

                    # update denominator
                    denominator *= ((Xfs[item['molecule']]
                                    * (P/P_ref)*(1)))**coefficient_

                # SECTION: loop over products
                for item in self.reaction_analysis[reaction]['products']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_]

                    # update numerator
                    numerator *= ((Xfs[item['molecule']]*(P/P_ref)*(phis[item['molecule']]))
                                  )**item['coefficient']

                # NOTE: update reaction term
                # reaction_terms[reaction] = (
                #     numerator/denominator) - equilibrium_data[reaction]

                reaction_terms[reaction] = numerator - \
                    denominator*Keq_dict[reaction]

            # SECTION: build objective functions
            obj = 0
            for i, reaction in enumerate(self.reaction_analysis):
                # NOTE: method 1
                # obj = abs(reaction_terms[reaction])
                # NOTE: method 2
                # obj += reaction_terms[reaction]**2
                obj += reaction_terms[reaction]**2

            # penalty obj
            # obj += 1e-3

            # sqrt
            obj = sqrt(obj)

            return obj
        except Exception as e:
            raise Exception(
                f"Error in ReactionOptimizer.equilibrium_reaction_objective_function(): {str(e)}") from e

    def reaction_equilibrium_equation(self,
                                      Xfs: Dict[str, float],
                                      P: float,
                                      T: float,
                                      Keq_dict: Dict[str, float],
                                      phase: Literal["liquid", "gas"] = "gas",
                                      gas_mixture: Literal["ideal",
                                                           "non-ideal"] = "ideal",
                                      liquid_mixture: Literal["ideal",
                                                              "non-ideal"] = "ideal",
                                      **kwargs):
        """
        Generate reaction equilibrium equations for gas and liquid phases.

        Parameters
        ----------
        Xfs : dict
            Mole fraction dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        P : float
            Pressure [Pa]
        T : float
            Temperature [K]
        Keq_dict : dict
            Reaction equilibrium constant dictionary calculated at T as: {key: value}, such as {R1: 0.1, R2: 0.2, ...}
        phase : str, optional
            Phase of calculation, by default "gas".
        gas_mixture : str, optional
            Mode for gas phase calculation, by default "ideal".
        liquid_mixture : str, optional
            Mode for liquid phase calculation, by default "ideal".
        kwargs : dict
            Additional parameters for non-ideal calculations.

        Notes
        -----
        The equilibrium constant at temperature T is calculated from the standard Gibbs free energy change:

        K_T = exp(-ΔG_T⁰ / RT) = ∏ [ (f̂_i / f_i⁰) ^ ν_i ]

        Where:
        - K_T    : Equilibrium constant at temperature T
        - ΔG_T⁰  : Standard Gibbs free energy change at temperature T (J/mol)
        - R      : Universal gas constant (8.314 J/mol·K)
        - T      : Temperature in Kelvin
        - f̂_i    : The mixture fugacity of component i.
        - f_i⁰   : The fugacity of pure component i in its reference state at temperature T.
        - ν_i    : Stoichiometric coefficient of component i (positive for products, negative for reactants)

        This expression connects thermodynamic data to fugacity-based equilibrium conditions.

        Gas Phase Calculation:
        - For gaseous compounds, the reference is the ideal gas at 1 bar, then f_i⁰ = P_ref = 1 bar.
        - f̂_i for non-ideal gases is calculated as f̂_i = φ_i * P * Y_i.
        - φ_i is the fugacity coefficient in the mixture.

        Liquid Phase Calculation:
        - f_i⁰ is standard-state fugacity (e.g., pure liquid at 1 bar)
        - f̂_i for liquid mixtures is defined as f̂_i = γ_i * X_i * f_i(pure,T).
        - f_i(pure,T) is the fugacity of pure component i at temperature T.
        - γ_i is the activity coefficient of component i in the liquid phase.
        - X_i is the mole fraction of component i in the liquid phase.

        Assumptions:
        - Regular liquid mixtures: Lewis-Randall
        - Ideal liquid solution: Raoult
        - Dilute or electrolyte systems: Henry
        - Poynting correction for liquid mixtures at high pressures.
        - General theoretical frameworks: Fugacity-based (chemical potential)
        """
        try:
            # SECTION: fugacity coefficient
            # fugacity coefficient
            fugacity_coeff = {}

            # check gas mixture
            if gas_mixture == "non-ideal":
                pass
            elif gas_mixture == "ideal":
                # set
                for i, key in enumerate(self.component_dict.keys()):
                    fugacity_coeff[key] = 1
            else:
                raise ValueError(
                    "Invalid gas mixture mode. Must be 'ideal' or 'non-ideal'.")

            # SECTION: activity coefficient
            # activity coefficient
            activity_coeff = {}

            # check liquid mixture
            if liquid_mixture == "non-ideal":
                pass
            elif liquid_mixture == "ideal":
                # set
                for i, key in enumerate(self.component_dict.keys()):
                    activity_coeff[key] = 1
            else:
                raise ValueError(
                    "Invalid liquid mixture mode. Must be 'ideal' or 'non-ideal'.")

            # SECTION: equilibrium equation
            # reaction term
            reaction_terms = {}

            # NOTE: loop over reactions
            for i, reaction in enumerate(self.reaction_analysis):
                # denominator
                denominator = 1
                # numerator
                numerator = 1

                # SECTION: loop over reactants
                for item in self.reaction_analysis[reaction]['reactants']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_]

                    # NOTE: cal
                    # check phase
                    if state_ == "gas":
                        # set
                        term_ = Xfs_ * P * fugacity_coeff[molecule_]
                    elif state_ == "liquid":
                        # set
                        # Lewis-Randall/Raoult
                        term_ = Xfs_ * activity_coeff[molecule_]
                    else:
                        raise ValueError(
                            "Invalid phase. Must be 'gas' or 'liquid'.")

                    # update denominator
                    denominator *= (term_)**coefficient_

                # SECTION: loop over products
                for item in self.reaction_analysis[reaction]['products']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_]

                    # NOTE: cal
                    # check phase
                    if state_ == "gas":
                        # set
                        term_ = Xfs_ * P * fugacity_coeff[molecule_]
                    elif state_ == "liquid":
                        # set
                        term_ = Xfs_ * activity_coeff[molecule_]
                    else:
                        raise ValueError(
                            "Invalid phase. Must be 'gas' or 'liquid'.")

                    # update numerator
                    numerator *= (term_)**coefficient_

                # NOTE: update reaction term
                # reaction_terms[reaction] = (
                #     numerator/denominator) - equilibrium_data[reaction]

                reaction_terms[reaction] = numerator - \
                    denominator*Keq_dict[reaction]
        except Exception as e:
            raise Exception(
                f"Error in ReactionOptimizer.equilibrium_reaction_objective_function(): {str(e)}") from e

    def constraint1(self, x, params):
        '''
        Sets a bound between 0 and 1 [x>>0 & x<<1] representing {-f(x) + 1 >= 0}

        Parameters
        ------------
        x : list
            loop values
        params : list
            parameters

        Returns
        -------
        obj : float
            objective function

        Notes
        -----
        1. component mole fraction should be less than 1
        '''
        # print(f"1: {x}")
        N0s, comp_list, component_dict, i = params

        # build EoR
        comp_value_list, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # constraint
        cons = -1*Xfs_vector[i] + 1

        return cons

    def constraint2(self, x, params):
        '''
        Sets a bound between 0 and 1 [x>>0 & x<<1] representing {f(x) >= 0}

        Parameters
        ------------
        x : list
        loop values
        params : list
        parameters

        Returns
        -------
        obj : float
        objective function

        Notes
        -----
        1. component mole fraction should be greater than 0
        2. define epsilon constraint 1e-5
        '''
        # print(f"2: {x}")
        N0s, comp_list, component_dict, i = params

        # build EoR
        comp_value_list, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # constraint
        cos = Xfs_vector[i] - 1e-5

        return cos

    def constraint3(self, x, params):
        '''
        Sets a bound for total mole f(x)>=0

        Parameters
        ------------
        x : list
        loop values
        params : list
        parameters

        Returns
        -------
        cons : float
        constraint

        Notes
        -----
        1. total mole of components should be greater than 0
        '''
        # print(f"3: {x}")
        N0s, comp_list, component_dict, i = params

        # build EoR
        _, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # cons
        cons = Nfs_vector[i]

        return cons

    def constraint4(self, x, params):
        '''
        Set a bound for reactants which are consumped ff(x)<=f0(x)

        Parameters
        ------------
        x : list
        loop values
        params : list
        parameters

        Returns
        -------
        cons : float
        constraint

        Notes
        -----
        1. final mole of reactants should be greater than 0
        '''
        # print(f"4: {x}")
        N0s, comp_list, component_dict, i = params

        # build EoR
        comp_value_list, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # cons
        cons = -1*Nfs_vector[i] + N0s_vector[i]

        return cons

    def constraint5(self, x, params):
        '''
        Sets a bound for sum of all mole fraction (x[i]<1)

        Parameters
        ------------
        x : list
        loop values
        params : list
        parameters

        Returns
        -------
        cons : float
        constraint

        Notes
        -----
        1. sum of all mole fraction should be 1
        '''
        # print(f"5: {x}")
        N0s, comp_list, component_dict = params

        # build EoR
        comp_value_list, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # cons
        cons = np.sum(Xfs_vector) - 1

        return cons

    def constraint6(self, x, params):
        '''
        Sets a bound for sum of all mole fraction (x[i] = 1) representing f(x)=1

        Parameters
        ------------
        x : list
        loop values
        params : list
        parameters

        Returns
        -------
        cons : float
        constraint

        Notes
        -----
        1. sum of all mole fraction should be 1
        '''
        # print(f"5: {x}")
        N0s, comp_list, component_dict = params

        # build EoR
        comp_value_list, comp_value_matrix = build_EoR(comp_list, x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # cons
        cons = np.sum(Xfs_vector) - 1

        return cons

    def opt_run(self, P, T, N0s, P_ref, EoR_init, component_dict, comp_list,
                analyze_overall_reactions_res, mode, reaction_res, Kas_T,
                fugacity_obj, phase, eos_model):

        # params
        params = (N0s, P_ref)

        # bounds
        bound0 = (0, 20)
        bounds = []
        for i in range(len(EoR_init)):
            bounds.append(bound0)

        # define constraint
        cons = []

        # inequality
        # looping through components
        for key, value in component_dict.items():
            # append constraint
            cons.append({'type': 'ineq', 'fun': constraint1, 'args': (
                (N0s, comp_list, component_dict, value),)})
            # append constraint
            cons.append({'type': 'ineq', 'fun': constraint2, 'args': (
                (N0s, comp_list, component_dict, value),)})
            # append constraint
            cons.append({'type': 'ineq', 'fun': constraint3, 'args': (
                (N0s, comp_list, component_dict, value),)})

            # check consumed
            if key in analyze_overall_reactions_res['consumed']:
                # append constraint
                cons.append({'type': 'ineq', 'fun': constraint4, 'args': (
                    (N0s, comp_list, component_dict, value),)})

        # append constraint
        cons.append({'type': 'ineq', 'fun': constraint5,
                    'args': ((N0s, comp_list, component_dict),)})

        # equality
        # constraint 5
        # cons.append({'type': 'ineq', 'fun': constraint5,'args':((N0s, comp_list,component_dict),)})

        # optimize
        opt_res = optimize.minimize(equilibrium_reaction_objective_function, EoR_init, method='SLSQP',
                                    args=(comp_list, reaction_res, Kas_T, params,
                                          P, T, component_dict, mode, eos_model, phase, fugacity_obj),
                                    bounds=bounds, constraints=cons, options={'disp': True, 'ftol': 1e-12, 'maxiter': 1000})

        # save
        return opt_res

    def process_optimization_results(self, res, input_data, comp_list, component_dict):
        '''
        Check optimization results

        Parameters
        ----------
        res : list
            optimization results
        initial_data : dict
            initial data
        comp_list : list
            component list, coefficient list
        component_dict : dict
            component dictionary

        Returns
        -------
        None
        '''
        # EoR
        EoR = res['x']

        # initial mole
        N0s = input_data['N0s']

        # unpack
        N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

        # update comp_list with EoR
        # build EoR => item[key] = value * EoR[i]
        _, comp_value_matrix = build_EoR(comp_list, EoR)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
            component_dict, N0s_vector, EoR_vector)

        # Nfs
        Nfs = {}
        for key, value in component_dict.items():
            Nfs[key] = Nfs_vector[value]

        return Xfs, Xfs_vector, Nf, Nfs_vector, Nfs
