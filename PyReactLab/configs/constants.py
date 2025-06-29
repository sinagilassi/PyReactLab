# CONSTANTS
# --------------
# import libs
import math

# NOTE: eos models
PENG_ROBINSON = "PR"
SOAVE_REDLICH_KWONG = "SRK"
REDLICH_KWONG = "RK"
VAN_DER_WAALS = "VDW"
# eos models
EOS_MODELS = [
    PENG_ROBINSON,
    SOAVE_REDLICH_KWONG,
    REDLICH_KWONG,
    VAN_DER_WAALS
]

# NOTE: assumptions
RAOULT_MODEL = 'raoult'
MODIFIED_RAOULT_MODEL = 'modified-raoult'

# NOTE: activity coefficient model
NRTL_ACTIVITY_MODEL = 'NRTL'
UNIQUAC_ACTIVITY_MODEL = 'UNIQUAC'
# activity models
ACTIVITY_MODELS = [
    NRTL_ACTIVITY_MODEL,
    UNIQUAC_ACTIVITY_MODEL
]

# NOTE: universal gas constant [J/mol.K]
R_CONST_J__molK = 8.314472

# NOTE: pi
PI_CONST = math.pi

# NOTE: STP condition
# pressure [Pa]
PRESSURE_STP_Pa = 101325
# temperature [K]
TEMPERATURE_STP_K = 273.15
# reference pressure [Pa]
PRESSURE_REF_Pa = 101325
# reference temperature [K]
TEMPERATURE_REF_K = 298.15

# SECTION: PyThermoDBLink/PyThermoDB
DATASOURCE = "datasource"
EQUATIONSOURCE = "equationsource"

# NOTE: set symbols
# enthalpy of formation ideal gas
EnFo_IG = "EnFo_IG"
# enthalpy of formation liquid
EnFo_LIQ = "EnFo_LIQ"
# gibbs free energy of formation ideal gas
GiEnFo_IG = "GiEnFo_IG"
# gibbs free energy of formation liquid
GiEnFo_LIQ = "GiEnFo_LIQ"

# SECTION: define symbols
# equilibrium constant
EQUILIBRIUM_CONSTANT_STD = 'equilibrium_constant_std'
EQUILIBRIUM_CONSTANT_STD_SYMBOL = 'K_eq_std'
EQUILIBRIUM_CONSTANT = 'equilibrium_constant'
EQUILIBRIUM_CONSTANT_SYMBOL = 'K_eq'
# Gibbs energy of reaction at 298.15 K
GIBBS_FREE_ENERGY_OF_REACTION_STD = 'gibbs_free_energy_of_reaction_std'
GIBBS_FREE_ENERGY_OF_REACTION_STD_SYMBOL = 'GiEn_rxn_std'
# enthalpy of reaction at 298.15 K
ENTHALPY_OF_REACTION_STD = 'enthalpy_of_reaction_std'
ENTHALPY_OF_REACTION_STD_SYMBOL = 'En_rxn_std'
# Gibbs energy of formation at 298.15 K
GIBBS_FREE_ENERGY_OF_FORMATION_STD = 'gibbs_free_energy_of_formation_std'
GIBBS_FREE_ENERGY_OF_FORMATION_STD_SYMBOL = 'GiEnFo'
# enthalpy of formation at 298.15 K
ENTHALPY_OF_FORMATION_STD = 'enthalpy_of_formation_std'
ENTHALPY_OF_FORMATION_STD_SYMBOL = 'EnFo'
# Gibbs energy of reaction at T
GIBBS_FREE_ENERGY_OF_REACTION_T = 'gibbs_free_energy_of_reaction_T'
GIBBS_FREE_ENERGY_OF_REACTION_T_SYMBOL = 'GiEn_rxn_T'
# enthalpy of reaction at T
ENTHALPY_OF_REACTION_T = 'enthalpy_of_reaction_T'
ENTHALPY_OF_REACTION_T_SYMBOL = 'En_rxn_T'
# gibbs energy of formation at T
GIBBS_FREE_ENERGY_OF_FORMATION_T = 'gibbs_free_energy_of_formation_T'
GIBBS_FREE_ENERGY_OF_FORMATION_T_SYMBOL = 'GiEnFo_T'
# enthalpy of formation at T
ENTHALPY_OF_FORMATION_T = 'enthalpy_of_formation_T'
ENTHALPY_OF_FORMATION_T_SYMBOL = 'EnFo_T'
# chemical potential in the mixture at T and P
CHEMICAL_POTENTIAL_MIXTURE_T_P = 'chemical_potential_mixture_T_P'
CHEMICAL_POTENTIAL_MIXTURE_T_P_SYMBOL = 'ChPo_mix_T_P'
# actual gibbs energy of reaction at T and P
ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION = 'actual_gibbs_free_energy_of_reaction'
ACTUAL_GIBBS_FREE_ENERGY_OF_REACTION_SYMBOL = 'GiEn_rxn_actual'
