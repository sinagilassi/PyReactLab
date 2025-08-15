# import packages/modules
from typing import (
    Dict,
    List,
    Literal,
    Optional,
    Any
)
import numpy as np
from math import sqrt, pow, exp, log
from scipy import optimize
from pyThermoModels import NRTL, UNIQUAC
from scipy.optimize import Bounds, NonlinearConstraint
# local
from ..configs import (
    R_CONST_J__molK,
    DATASOURCE,
    EQUATIONSOURCE,
    PRESSURE_REF_Pa,
    TEMPERATURE_REF_K,
    EOS_MODELS,
    ACTIVITY_MODELS
)


class ReactionOptimizer:
    '''Reaction Optimizer class'''

    # SECTION: variables
    _system_inputs = None
    # universal gas constant [J/mol.K]
    R = R_CONST_J__molK
    # temperature [K]
    T_Ref_K = TEMPERATURE_REF_K
    # pressure [Pa]
    P_Ref_Pa = PRESSURE_REF_Pa
    # pressure [bar]
    P_Ref_bar = 1

    # NOTE: models
    # eos
    _eos = None
    # eos model
    _eos_model = None
    # activity model
    _activity_model = None
    # activity
    _activity = None
    # activity inputs
    _activity_inputs = None

    # NOTE: systems
    # gas phase
    _gas_mixture = 'ideal'
    # solution
    _solution = 'ideal'

    def __init__(
        self,
        datasource: Dict[str, Any],
        equationsource: Dict[str, Any],
        component_dict: Dict[str, float | int],
        comp_list: List[Dict[str, float | int]],
        stoichiometric_coeff: np.ndarray,
        reaction_analysis: Dict,
        overall_reaction_analysis: Dict,
        **kwargs
    ):
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
        # set stoichiometric coefficient
        self.stoichiometric_coeff = stoichiometric_coeff
        # set reaction analysis
        self.reaction_analysis = reaction_analysis
        # set overall reaction analysis
        self.overall_reaction_analysis = overall_reaction_analysis

        # NOTE: set
        # component list
        self.component_list = list(component_dict.keys())

        # NOTE: kwargs
        self.threshold = kwargs.get('threshold', 1e-8)

    @property
    def eos(self):
        '''Get eos class'''
        # check
        return self._eos

    @eos.setter
    def eos(self, value: Any):
        '''Set eos class'''
        # set
        self._eos = value

    @property
    def eos_model(self):
        '''Get eos model'''
        return self._eos_model

    @eos_model.setter
    def eos_model(self, value: str):
        '''Set eos model'''
        # set
        self._eos_model = value

    @property
    def activity_model(self):
        '''Get activity model'''
        return self._activity_model

    @activity_model.setter
    def activity_model(self, value: str):
        '''Set activity model'''
        # set
        self._activity_model = value

    @property
    def activity(self):
        '''activity model (NRTL, UNIQUAC)'''
        return self._activity

    @activity.setter
    def activity(self, value: Any):
        '''Set activity model'''
        # set
        self._activity = value

    @property
    def activity_inputs(self):
        '''Get activity inputs'''
        return self._activity_inputs

    @activity_inputs.setter
    def activity_inputs(self, value: Dict[str, Any]):
        '''Set activity inputs'''
        # set
        self._activity_inputs = value

    @property
    def gas_mixture(self):
        '''Get gas mixture'''
        return self._gas_mixture

    @gas_mixture.setter
    def gas_mixture(self, value: str):
        '''Set gas mixture'''
        # set
        self._gas_mixture = value

    @property
    def solution(self):
        '''Get solution'''
        return self._solution

    @solution.setter
    def solution(self, value: str):
        '''Set solution'''
        # set
        self._solution = value

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

            # ! looping through reactions
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

    def build_final_X(
        self,
        N0s_vector: np.ndarray,
        EoR_vector: np.ndarray
    ):
        '''
        build final X

        Parameters
        ----------
        component_dict : dict
            component_dict
        N0s_vector : numpy.ndarray
            N0s_vector
        EoR_vector : numpy.ndarray
            EoR_vector

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

    def obj_fn(
        self,
        x,
        N0s: Dict[str, float],
        P: float,
        T: float,
        equilibrium_constants: Dict[str, float],
        method: Literal['minimize', 'least_squares'],
    ):
        '''
        Objective function for optimization

        Parameters
        ----------
        x : list
            loop values
        N0s : dict
            initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
            parameters
        P : float
            pressure [bar]
        T : float
            temperature [K]
        equilibrium_constants : dict
            reaction equilibrium constant dictionary calculated at T as: {key: value}, such as {R1: 0.1, R2: 0.2, ...}
        method : str
            optimization method.

        Returns
        -------
        obj : float | list
            objective function, if method is 'minimize' return float, if method is 'least_squares' return list
        '''
        try:
            # NOTE: default values
            # pressure [bar]
            # self.P_Ref_bar

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

            # SECTION: equilibrium equation
            # reaction term
            reaction_terms = self.reaction_equilibrium_equation(
                Xfs=Xfs,
                P=P,
                T=T,
                equilibrium_constants=equilibrium_constants,
            )

            # SECTION: build objective functions
            if method == 'minimize':
                obj = 0
                for i, reaction in enumerate(self.reaction_analysis):
                    # NOTE: method 1
                    # obj = abs(reaction_terms[reaction])
                    # NOTE: method 2
                    # obj += reaction_terms[reaction]**2
                    # NOTE: set
                    obj += reaction_terms[reaction]**2

                # NOTE: define objective function
                # penalty obj
                # obj += 1e-3

                # sqrt
                # obj = sqrt(obj)
            elif method == 'least_squares':
                # NOTE: set
                obj = []
                for i, reaction in enumerate(self.reaction_analysis):
                    # NOTE: method 1
                    # obj.append(reaction_terms[reaction])
                    # NOTE: method 2
                    # obj.append(reaction_terms[reaction]**2)
                    # NOTE: set
                    obj.append(reaction_terms[reaction])

                # NOTE: mole fraction constraints
                Xfs_sum = np.sum(Xfs_vector)
                obj.append(Xfs_sum - 1)
            else:
                raise ValueError(
                    f"Invalid optimization method: {method}. Must be 'minimize' or 'least_squares'.")

            return obj
        except Exception as e:
            raise Exception(
                f"Error in objective function: {str(e)}") from e

    def reaction_equilibrium_equation(
        self,
        Xfs: Dict[str, float],
        P: float,
        T: float,
        equilibrium_constants: Dict[str, Dict[str, Any]],
        phase: Literal[
            "liquid", "gas"
        ] = "gas",
        **kwargs
    ):
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
        equilibrium_constants : dict
            Reaction equilibrium constant dictionary calculated at T as: such as {R1: dict, R2: dict, ...},
            the dict value consists of `value`, `symbol`, `unit`, `temperature`, `reaction`, `method`
        phase : str, optional
            Phase of calculation, by default "gas".
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
            # NOTE: default values
            # pressure [bar]
            # self.P_Ref_bar
            # temperature [K]
            # self.T_Ref_K

            # SECTION: fugacity coefficient
            # fugacity coefficient
            fugacity_coeff = {}

            # check gas mixture
            if self.gas_mixture.lower() == "non-ideal":
                # NOTE: model input
                model_input = {
                    "feed-specification": Xfs,
                    "pressure": [P, 'bar'],
                    "temperature": [T, 'K'],
                }

                # NOTE: eos model
                eos_model = self.eos_model

                # Validate the model name
                if eos_model not in EOS_MODELS:
                    raise ValueError(
                        f"Invalid EOS model: {eos_model}. Must be {EOS_MODELS}.")

                # NOTE: calculate fugacity
                res_ = self._cal_fugacity_coefficient_gaseous_mixture(
                    model_name=eos_model,  # type: ignore
                    model_input=model_input)
                # update
                fugacity_coeff = {**res_}

            elif self.gas_mixture.lower() == "ideal":
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
            if self.solution.lower() == "non-ideal":

                # prepare model input
                # Ensure activity_inputs is a dictionary before unpacking
                activity_inputs = self.activity_inputs or {}
                # set
                model_inputs = {
                    "mole_fraction": Xfs, **activity_inputs
                }

                # NOTE: check model name
                if self.activity_model == 'NRTL':
                    # exec
                    activity_coeff = self._cal_activity_coefficient_solution(
                        model_name='NRTL', model_input=model_inputs)
                elif self.activity_model == 'UNIQUAC':
                    # exec
                    activity_coeff = self._cal_activity_coefficient_solution(
                        model_name='UNIQUAC', model_input=model_inputs)
                else:
                    raise ValueError(
                        f"Invalid activity model: {self.activity_model}. Must be {ACTIVITY_MODELS}.")
            elif self.solution.lower() == "ideal":
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

                # state count
                state_count_ = self.reaction_analysis[reaction]['state_count']

                # SECTION: loop over reactants
                for item in self.reaction_analysis[reaction]['reactants']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    molecule_state_ = item['molecule_state']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_state_]

                    # NOTE: cal
                    # check phase
                    if state_ == "g":
                        # set
                        term_ = Xfs_ * (P / self.P_Ref_bar) * \
                            fugacity_coeff[molecule_state_]
                    elif state_ == "l":
                        # check solution
                        if self.solution.lower() == "non-ideal":
                            # Lewis-Randall/Raoult
                            term_ = Xfs_ * activity_coeff[molecule_state_]
                        elif self.solution.lower() == "ideal":
                            # check
                            if state_count_['l'] == 1:
                                # pure liquid
                                term_ = 1
                        else:
                            raise ValueError(
                                "Invalid liquid mixture mode. Must be 'ideal' or 'non-ideal'.")

                    elif state_ == "s":
                        # set
                        # solid always has the activity of 1
                        term_ = 1
                    else:
                        raise ValueError(
                            "Invalid phase. Must be 'gas' or 'liquid'.")

                    # update denominator
                    denominator *= (term_)**coefficient_

                # SECTION: loop over products
                for item in self.reaction_analysis[reaction]['products']:
                    # NOTE: item info
                    molecule_ = item['molecule']
                    molecule_state_ = item['molecule_state']
                    coefficient_ = item['coefficient']
                    state_ = item['state']

                    # final mole fraction
                    Xfs_ = Xfs[molecule_state_]

                    # NOTE: cal
                    # check phase
                    if state_ == "g":
                        # set
                        term_ = Xfs_ * (P / self.P_Ref_bar) * \
                            fugacity_coeff[molecule_state_]
                    elif state_ == "l":
                        # check solution
                        if self.solution.lower() == "non-ideal":
                            # Lewis-Randall/Raoult
                            term_ = Xfs_ * activity_coeff[molecule_state_]
                        elif self.solution.lower() == "ideal":
                            # check
                            if state_count_['l'] == 1:
                                # pure liquid
                                term_ = 1
                        else:
                            raise ValueError(
                                "Invalid liquid mixture mode. Must be 'ideal' or 'non-ideal'.")
                    elif state_ == "s":
                        # set
                        # solid always has the activity of 1
                        term_ = 1
                    else:
                        raise ValueError(
                            "Invalid phase. Must be 'gas' or 'liquid'.")

                    # update numerator
                    numerator *= (term_)**coefficient_

                # NOTE: update reaction term
                # equilibrium data
                Keq_ = equilibrium_constants[reaction]['value']

                # ? method 1
                # reaction_terms[reaction] = (numerator/denominator) - Keq_
                # ? method 2
                # reaction_terms[reaction] = numerator - denominator*Keq_
                # ? method 3
                Q_i = numerator / denominator
                Q_i_safe = np.clip(Q_i, 1e-12, None)
                reaction_terms[reaction] = log(Q_i_safe) - log(Keq_)

            # return
            return reaction_terms
        except Exception as e:
            raise Exception(
                f"Error in building the reaction equilibrium equation for {phase}: {str(e)}") from e

    def _cal_fugacity_coefficient_gaseous_mixture(
        self,
        model_name: Literal[
            "SRK", "PR", "RK"
        ],
        model_input: Dict[str, Any]
    ):
        """
        Calculate the fugacity coefficient of gaseous mixture using the specified EOS model.

        Parameters
        ----------
        model_name : str
            The EOS model to use for calculation. Options: "SRK", "PR", "RK".
        model_input : dict
            The input data for the EOS model.
            - feed-specification: Dictionary of component mole fractions.
            - pressure: pressure value in any unit.
            - temperature: temperature value in any unit.

        Returns
        -------
        dict
            A dictionary containing the fugacity coefficients for each component in the gaseous mixture.
        """
        try:
            # SECTION: model source
            model_source = {
                "datasource": self.datasource,
                "equationsource": self.equationsource
            }

            # NOTE: calculate fugacity
            eos = self.eos

            # check
            if eos is None or eos == 'None':
                raise ValueError(
                    f"Invalid EOS model. Must be {EOS_MODELS}.")

            res = eos.cal_fugacity_mixture(
                model_name=model_name,
                model_input=model_input,
                model_source=model_source
            )

            # NOTE: extract fugacity coefficient for each component
            res_1 = res['vapor']

            # res
            phi_comp = {}

            # looping through each component
            for i, key in enumerate(res_1.keys()):
                # set
                phi_comp[key] = res_1[key]['fugacity_coefficient']['value']

            return phi_comp
        except Exception as e:
            raise Exception(
                f"Error in calculating the fugacity coefficient for the gaseous mixture: {str(e)}") from e

    def _cal_activity_coefficient_solution(
        self,
        model_name: Literal[
            'NRTL', 'UNIQUAC'
        ],
        model_input: Dict[str, Any],
    ):
        """
        Calculate the activity coefficient of solution using the specified model.

        Parameters
        ----------
        model_name : str
            The model to use for calculation. Options: "NRTL", "UNIQUAC".
        model_input : dict
            The input data for the model.
            - mole_fraction: Dictionary of component mole fractions.
            - tau_ij: Dictionary of interaction parameters.
            - alpha_ij: Dictionary of interaction parameters.
            - r_i: relative van der Waals volume of component i
            - q_i: relative surface area of component i
        model_source : dict
            The source of the model data.
            - datasource: Data source for the model.
            - equationsource: Equation source for the model.

        Returns
        -------
        dict
            A dictionary containing the fugacity coefficients for each component in the liquid mixture.
        """
        try:
            # SECTION: model source
            # model_source = {
            #     "datasource": self.datasource,
            #     "equationsource": self.equationsource
            # }
            # check activity
            if self.activity is None or self.activity == 'None':
                raise ValueError(
                    f"Invalid activity model. Must be {ACTIVITY_MODELS}.")

            # NOTE: model name
            if model_name == "NRTL":
                # calculate fugacity
                res, others = self.activity.cal(model_input=model_input)
            elif model_name == "UNIQUAC":
                # calculate fugacity
                res, others = self.activity.cal(model_input=model_input)
            else:
                raise ValueError(
                    f"Invalid model name: {model_name}. Must be 'NRTL' or 'UNIQUAC'.")

            # NOTE: extract fugacity coefficient for each component
            return others['AcCo_i_comp']
        except Exception as e:
            raise Exception(
                f"Error in calculating the fugacity coefficient for the liquid mixture: {str(e)}") from e

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
        comp_value_list, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

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
        comp_value_list, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

        # constraint
        # ! set epsilon constraint
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
        _, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

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
        comp_value_list, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

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
        comp_value_list, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

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
        comp_value_list, comp_value_matrix = self.build_EoR(x)

        # extent sun [mol]
        EoR_vector = np.sum(comp_value_matrix, axis=0)

        # unpack N0s
        N0s_list, N0s_vector, N0f = self.unpack_X(N0s)

        # build final X
        Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
            N0s_vector, EoR_vector)

        # cons
        cons = np.sum(Xfs_vector) - 1

        return cons

    def opt_run(
        self,
        initial_mole: Dict[str, float | int],
        initial_mole_fraction: Dict[str, float | int],
        pressure: float,
        temperature: float,
        equilibrium_constants: Dict[str, float],
        reaction_numbers: int,
        method: Literal['minimize', 'least_squares'] = 'minimize',
        **kwargs
    ):
        """
        Start optimization process

        Parameters
        ----------
        initial_mole : dict
            Initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        initial_mole_fraction : dict
            Initial mole fraction dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        pressure : float
            Pressure [bar]
        temperature : float
            Temperature [K]
        equilibrium_constants : dict
            Reaction equilibrium constant dictionary calculated at T as: {key: value}, such as {R1: 0.1, R2: 0.2, ...}
        reaction_numbers : int
            Number of reactions
        method : str, optional
            Optimization method, by default 'minimize'.
        kwargs : dict
            Additional parameters for optimization.
            - least_square_algorithm: 'trf', 'dogbox', 'lm'
            - minimize_algorithm: 'SLSQP', 'trust-constr', 'COBYLA'
            - scale: Scale factor for the upper bound.
            - bound_scale: Scale factor for the upper bound.
            - fallback: Fallback value for upper bound if no reactants are present.

        Returns
        -------
        opt_res : OptimizeResult
            Optimization result object containing the optimized values and other information.
        """
        try:
            # NOTE: default values
            # least_square_algorithm
            least_square_algorithm = kwargs.get(
                'least_square_algorithm', 'trf')
            # minimize_algorithm
            minimize_algorithm = kwargs.get('minimize_algorithm', 'SLSQP')

            # ? check if the algorithm is valid
            if minimize_algorithm not in ['SLSQP', 'trust-constr', 'COBYLA']:
                raise ValueError(
                    f"Invalid minimize algorithm: {minimize_algorithm}. Must be 'SLSQP', 'trust-constr', or 'COBYLA'.")

            # ? check if the algorithm is valid
            if least_square_algorithm not in ['trf', 'dogbox', 'lm']:
                raise ValueError(
                    f"Invalid least square algorithm: {least_square_algorithm}. Must be 'trf', 'dogbox', or 'lm'.")

            # NOTE: set
            # initial guess for extent of reaction
            # default value
            EOR0_val = kwargs.get('EOR0_val', 0.5)

            # set
            EOR0 = np.random.uniform(0, EOR0_val, reaction_numbers)

            # NOTE: bounds
            bound0 = (0, 20)
            bounds = []
            for i in range(len(EOR0)):
                bounds.append(bound0)

            # EOR02
            # initial mole
            initial_mole_ = np.array(list(initial_mole.values()))

            # extent of reaction (EoR)
            EOR0_Bounds, EOR0_bounds, EoR_initial = self.compute_bounds(
                nu=self.stoichiometric_coeff,
                n0=initial_mole_,
                **kwargs,
            )

            # NOTE: constraints
            # constraints 1
            cons1 = self.constraints_collection_1(initial_mole)
            # constraints 2
            cons2 = self.constraints_collection_2(
                nu=self.stoichiometric_coeff,
                n0=initial_mole_,
                include_mole_fraction_constraint=True
            )

            # set constraints
            cons = []

            # SECTION: optimize
            if method == 'minimize':
                opt_res = optimize.minimize(
                    fun=self.obj_fn,
                    x0=EOR0,
                    args=(
                        initial_mole,
                        pressure,
                        temperature,
                        equilibrium_constants,
                        method
                    ),
                    method=minimize_algorithm,
                    bounds=EOR0_bounds,
                    constraints=cons1,
                    options={
                        'disp': False,
                        'ftol': 1e-12,
                        'maxiter': 1000
                    }
                )
            elif method == 'least_squares':
                opt_res = optimize.least_squares(
                    fun=self.obj_fn,
                    x0=EOR0,
                    args=(
                        initial_mole,
                        pressure,
                        temperature,
                        equilibrium_constants,
                        method
                    ),
                    bounds=EOR0_Bounds,
                    method=least_square_algorithm,
                )
            else:
                raise ValueError(
                    f"Invalid optimization method: {method}. Must be 'minimize' or 'least_squares'.")

            # save
            return opt_res
        except Exception as e:
            raise Exception(
                f"Error in the optimization process: {str(e)}") from e

    def constraints_collection_1(self, initial_mole: Dict[str, float | int]):
        """
        Generate a collection of constraints for optimization.

        Parameters
        ----------
        initial_mole : Dict[str, float | int]
            Initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        """
        try:
            # NOTE: define constraint
            cons = []

            # inequality
            # looping through components
            for key, value in self.component_dict.items():
                # append constraint
                cons.append({
                    'type': 'ineq',
                    'fun': self.constraint1,
                    'args': (
                        (initial_mole, self.comp_list, self.component_dict, value),)
                })
                # append constraint
                cons.append({
                    'type': 'ineq',
                    'fun': self.constraint2,
                    'args': (
                        (initial_mole, self.comp_list, self.component_dict, value),)
                })
                # append constraint
                cons.append({
                    'type': 'ineq',
                    'fun': self.constraint3,
                    'args': (
                        (initial_mole, self.comp_list, self.component_dict, value),)
                })

                # check consumed
                if key in self.overall_reaction_analysis['consumed']:
                    # append constraint
                    cons.append({
                        'type': 'ineq',
                        'fun': self.constraint4,
                        'args': (
                            (initial_mole, self.comp_list, self.component_dict, value),)
                    })

            # append constraint
            cons.append({
                'type': 'ineq',
                'fun': self.constraint5,
                        'args': ((initial_mole, self.comp_list, self.component_dict),)
            })

            # equality
            # constraint 5
            # cons.append({
            # 'type': 'ineq',
            # 'fun': constraint5,
            # 'args':((N0s, comp_list,component_dict),)
            # })

            # return
            return cons
        except Exception as e:
            raise Exception(
                f"Error in generating constraints collection: {str(e)}") from e

    def constraints_collection_2(
        self,
        nu,
        n0,
        include_mole_fraction_constraint=True
    ):
        """
        Generate a collection of constraints for optimization.

        Parameters
        ----------
        initial_mole : Dict[str, float | int]
            Initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        """
        try:
            n_species, n_reactions = nu.shape
            constraints = []

            # 1. Species non-negativity constraints
            for i in range(n_species):
                def make_fn(i):
                    return lambda ksi: n0[i] + np.dot(nu[i, :], ksi)
                constraints.append(NonlinearConstraint(make_fn(i), 0, np.inf))

            # 2. Mole fraction sum constraint (optional)
            if include_mole_fraction_constraint:
                def mole_fraction_sum(ksi):
                    n = n0 + nu @ ksi
                    total = np.sum(n)
                    # avoid division by zero if total is very small
                    if total <= 1e-12:
                        return 0.0
                    y = n / total
                    return np.sum(y)

                constraints.append(NonlinearConstraint(
                    mole_fraction_sum, 1.0, 1.0))

            return constraints
        except Exception as e:
            raise Exception(
                f"Error in generating constraints collection: {str(e)}") from e

    def equilibrium_results(
        self,
        initial_mole: Dict[str, float | int],
        EoR: list
    ):
        '''
        Calculate the equilibrium results based on the initial mole and extent of reaction.

        Parameters
        ----------
        initial_mole : dict
            Initial mole dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        EoR : list
            Extent of reaction list [EoR1, EoR2, EoR3, ...]

        Returns
        -------
        Xfs : dict
            Final mole fraction dictionary {key: value}, such as {CO2: 0.1, H2: 0.1, CO: 0.2, H2O: 0.3, CH3OH: 0.4}
        Xfs_vector : np.ndarray
            Final mole fraction vector [Xfs1, Xfs2, Xfs3, ...]
        Nf : float
            Final moles of the system
        Nfs_vector : np.ndarray
            Final moles vector [Nf1, Nf2, Nf3, ...]
        Nfs : dict
            Final moles dictionary {key: value}, such as {CO2: 0, H2: 1, CO: 2, H2O: 3, CH3OH: 4}
        '''
        try:
            # unpack
            N0s_list, N0s_vector, N0f = self.unpack_X(initial_mole)

            # update comp_list with EoR
            # build EoR => item[key] = value * EoR[i]
            _, comp_value_matrix = self.build_EoR(EoR)

            # extent sun [mol]
            EoR_vector = np.sum(comp_value_matrix, axis=0)

            # build final X
            Xfs, Xfs_vector, Nf, Nfs_vector = self.build_final_X(
                N0s_vector, EoR_vector
            )

            # Nfs
            Nfs = {}
            for key, value in self.component_dict.items():
                Nfs[key] = Nfs_vector[value]

            # set float
            Nfs = {k: float(v) for k, v in Nfs.items()}
            Xfs = {k: float(v) for k, v in Xfs.items()}

            return Xfs, Xfs_vector, Nf, Nfs_vector, Nfs
        except Exception as e:
            raise Exception(
                f"Error in processing optimization results: {str(e)}") from e

    def compute_bounds(
        self,
        nu: np.ndarray,
        n0: np.ndarray,
        **kwargs
    ):
        """
        Computes lower and upper bounds for each reaction extent ξ_j.

        Parameters
        ----------
        nu: np.ndarray
            The stoichiometric matrix (n_species x n_reactions)
        N0s: np.ndarray
            Initial moles of each species (length: n_species)
        kwargs: dict
            Additional parameters for bounds calculation.
            - scale: float
                Scale factor for the upper bound.
            - bound_scale: float
                Scale factor for the upper bound.
            - fallback: float
                Fallback value for upper bound if no reactants are present.

        Returns
        -------
        Bounds
            scipy.optimize.Bounds object for use in minimize/least_squares.
        """
        try:
            # SECTION: default values
            # scale
            scale = kwargs.get('scale', 10.0)
            # bound_scale
            bound_scale = kwargs.get('bound_scale', 0.5)
            # fallback
            fallback = kwargs.get('fallback', None)

            # SECTION: Check if nu is a 2D array
            if nu.ndim != 2:
                raise ValueError("nu must be a 2D array.")
            # Check if n0 is a 1D array
            if n0.ndim != 1:
                raise ValueError("n0 must be a 1D array.")

            # NOTE: Check if nu and n0 have compatible dimensions
            if nu.shape[0] != n0.shape[0]:
                raise ValueError("nu and n0 must have compatible dimensions.")

            # NOTE: set fallback
            if fallback is None:
                fallback = max(n0.sum(), 1.0) * scale  # 10× total initial mol

            # SECTION: Compute bounds
            n_species, n_reactions = nu.shape
            lb = np.zeros(n_reactions)
            ub = np.full(n_reactions, np.inf)

            # Loop over reactions
            for j in range(n_reactions):
                for i in range(n_species):
                    v = nu[i, j]
                    if v < 0:
                        if n0[i] > 0:
                            max_mu = n0[i] / abs(v)
                            ub[j] = min(ub[j], max_mu)
                        else:
                            # No initial amount of a required reactant → can't proceed
                            ub[j] = min(ub[j], fallback)

            # NOTE: convert to bounds list
            bounds = []
            for i in range(len(lb)):
                # add slack to the upper bound
                ub[i] = ub[i]
                bounds.append((lb[i], ub[i]))

            # NOTE: initial guess
            EOR0 = []
            for i in range(len(ub)):
                # check if upper bound is finite
                if ub[i] == scale:
                    # set
                    _val = lb[i] + bound_scale * (ub[i] - lb[i])
                    EOR0.append(_val)
                else:
                    # set
                    EOR0.append(np.random.uniform(lb[i], ub[i]))

            return Bounds(lb, ub), bounds, EOR0
        except Exception as e:
            raise Exception(
                f"Error in computing bounds: {str(e)}") from e
