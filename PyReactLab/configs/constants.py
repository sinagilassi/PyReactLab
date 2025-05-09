# CONSTANTS
# --------------
# import libs
import math

# NOTE: eos models
PENG_ROBINSON = "PR"
SOAVE_REDLICH_KWONG = "SRK"
REDLICH_KWONG = "RK"
VAN_DER_WAALS = "VDW"

# NOTE: assumptions
RAOULT_MODEL = 'raoult'
MODIFIED_RAOULT_MODEL = 'modified-raoult'

# NOTE: activity coefficient model
NRTL_ACTIVITY_MODEL = 'NRTL'
UNIQUAC_ACTIVITY_MODEL = 'UNIQUAC'

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
