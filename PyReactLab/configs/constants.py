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
