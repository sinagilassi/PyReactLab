# import libs
from typing import Dict, Any, Literal, Optional


def cal_conversion(N0s, Nfs, component_item):
    '''
    Calculate conversion

    Parameters
    ----------
    input_data : dict
        input data
    Nfs : dict
        final mole
    component_item : str
        component item

    Returns
    -------
    conversion_res : dict
        conversion result
    '''
    # conversion
    # initial mole
    _N0i = N0s[component_item]
    # final mole
    _Nfi = Nfs[component_item]
    # conversion
    conversion = ((_N0i - _Nfi)/_N0i)*100

    return conversion


def cal_selectivity(input_data, Nfs, component_item, component_ref, carbon_count):
    '''
    Calculate selectivity

    Parameters
    ----------
    input_data : dict
        input data
    Nfs : dict
        final mole
    component_item : str
        component item
    component_ref : str
        component reference

    Returns
    -------
    selectivity : float
        selectivity result
    '''

    # selectivity to [mol]
    N0i_ref = input_data['N0s'][component_ref]
    Nfi_ref = Nfs[component_ref]

    # initial mole
    _N0i = input_data['N0s'][component_item]
    # final mole
    _Nfi = Nfs[component_item]
    # selectivity
    selectivity = (carbon_count*_Nfi - carbon_count*_N0i)/(N0i_ref - Nfi_ref)
    selectivity = selectivity*100

    # res
    return selectivity


def cal_yield(input_data, Nfs, component_item, component_ref, carbon_count=1):
    '''
    Calculate yield

    Parameters
    ----------
    input_data : dict
        input data
    Nfs : dict
        final mole
    component_item : str
        component item
    component_ref : str
        component reference
    carbon_count : int
        carbon count

    Returns
    -------
    yield_res : float
        yield result
    '''
    # reference
    N0i_ref = input_data['N0s'][component_ref]

    # component
    _N0i = input_data['N0s'][component_item]
    _Nfi = Nfs[component_item]

    # yield
    yield_res = ((_Nfi-_N0i)/N0i_ref)*100

    return yield_res
