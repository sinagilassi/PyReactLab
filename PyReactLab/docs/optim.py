# import packages/modules
import numpy as np
from math import sqrt, pow, exp
from scipy import optimize


def unpack_X(mol_data_pack, component_dict):
    '''
    convert X to dict {key: value} and list [value]

    Parameters
    ----------
    mol_data_pack : dict
        # mol_data_pack

    Returns
    -------
    N0s_list : list
        N0s_list
    N0s_vector : numpy.ndarray
        N0s_vector
    Nf : float
        Nf
    '''
    # size
    size = len(mol_data_pack)
    # lists
    N0s_list = []
    N0s_vector = np.zeros(size)

    # looping through
    for i, (com_key, com_value) in enumerate(component_dict.items()):
        N0s_list.append(mol_data_pack[com_key])
        N0s_vector[i] = mol_data_pack[com_key]

    # final
    Nf = np.sum(N0s_vector)

    return N0s_list, N0s_vector, Nf


def build_EoR(comp_list, EoR):
    '''
    build EoR

    Parameters
    ----------
    comp_list : list
        comp_list
    EoR : list
        EoR

    Returns
    -------
    comp_value_list : list
        comp_value_list
    comp_value_matrix : numpy.ndarray
        comp_value_matrix
    '''
    # comp value list
    comp_value_list = []
    # new
    comp_list_updated = []
    # copy data
    for i, item in enumerate(comp_list):
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


def build_final_X(component_dict, N0s_vector, EoR_vector):
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
    # final mole of components [mole]
    Nfs_vector = N0s_vector + EoR_vector

    # final mole fraction
    Xfs_vector = Nfs_vector/np.sum(Nfs_vector)

    # final mole
    Nf = np.sum(Nfs_vector)

    # convert to Xfs
    Xfs = {}
    for i, key in enumerate(component_dict.keys()):
        Xfs[key] = Xfs_vector[i]

    return Xfs, Xfs_vector, Nf, Nfs_vector


def equilibrium_reaction_objective_function(x, comp_list, reactions, equilibrium_data, params,
                                            P, T, component_dict, mode="ideal", eos_model='SRK', phase='VAPOR',
                                            fugacity_obj=None, threshold=1e-8):
    '''
    Equilibrium Reaction Analysis

    Parameters
    ----------
    x : list
        loop values
    comp_list : list
        component list
    reactions : list
        reaction list
    equilibrium_data : dict
        equilibrium data
    params : list
        parameters
    P : float
        pressure
    T : float
        temperature
    mode : str, optional
        mode of analysis, by default "ideal"

    Returns
    -------
    obj : float
        objective function
    '''
    # extent of reaction
    EoR = x

    # extent of reaction
    EoR = x
    # define threshold
    for i, item in enumerate(EoR):
        EoR[i] = max(threshold, item)

    # params
    N0s, P_ref = params

    # unpack N0s
    N0s_list, N0s_vector, N0f = unpack_X(N0s, component_dict)

    # update comp_list with EoR
    # build EoR => item[key] = value * EoR[i]
    _, comp_value_matrix = build_EoR(comp_list, EoR)

    # extent sun [mol]
    EoR_vector = np.sum(comp_value_matrix, axis=0)

    # build final X
    Xfs, Xfs_vector, Nf, Nfs_vector = build_final_X(
        component_dict, N0s_vector, EoR_vector)

    # calculate fugacity coefficient
    # model input
    # component list
    comp_list = list(component_dict.keys())

    # feed spec
    feed_spec = {}
    for comp in comp_list:
        feed_spec[comp] = Xfs_vector[comp_list.index(comp)]

    # model input
    model_input = {
        "eos-model": eos_model,
        "phase": phase,
        "feed-spec": feed_spec,
        "operating-conditions": {
            "pressure": [P, 'Pa'],
            "temperature": [T, 'K'],
        },
    }

    # calculate fugacity coefficient
    phis = []
    # check
    if mode == "non-ideal":
        # fugacity coefficient [-]
        fugacity_res = fugacity_obj.cal_fugacity_coefficient(model_input)
        # phii
        _, _, _, phi_pack = fugacity_res

        # phis
        phis = {}
        for i, key in enumerate(component_dict.keys()):
            phis[key] = phi_pack['VAPOR'][key]['phi']

    elif mode == 'ideal':
        phis = {}
        for i, key in enumerate(component_dict.keys()):
            phis[key] = 1
    else:
        raise ValueError("Invalid mode. Must be 'ideal' or 'non-ideal'.")

    # reaction term
    reaction_terms = {}

    for i, reaction in enumerate(reactions):
        # denominator
        denominator = 1
        # numerator
        numerator = 1
        # loop over reactants
        for item in reactions[reaction]['reactants']:
            # update denominator
            denominator *= ((Xfs[item['molecule']]*(P/P_ref)*(phis[item['molecule']]))
                            )**item['coefficient']
        # loop over products
        for item in reactions[reaction]['products']:
            # update numerator
            numerator *= ((Xfs[item['molecule']]*(P/P_ref)*(phis[item['molecule']]))
                          )**item['coefficient']
        # update reaction term
        # reaction_terms[reaction] = (
        #     numerator/denominator) - equilibrium_data[reaction]
        reaction_terms[reaction] = numerator - \
            denominator*equilibrium_data[reaction]

    # build objective functions
    obj = 0
    for i, reaction in enumerate(reactions):
        # method 1
        # obj = abs(reaction_terms[reaction])
        # method 2
        # obj += reaction_terms[reaction]**2
        obj += reaction_terms[reaction]**2

    # penalty obj
    # obj += 1e-3

    # sqrt
    obj = sqrt(obj)

    return obj


def constraint1(x, params):
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


def constraint2(x, params):
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


def constraint3(x, params):
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


def constraint4(x, params):
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


def constraint5(x, params):
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


def constraint6(x, params):
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


def opt_run(P, T, N0s, P_ref, EoR_init, component_dict, comp_list,
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


def process_optimization_results(res, input_data, comp_list, component_dict):
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
