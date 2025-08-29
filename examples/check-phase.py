# import libs
import PyReactLab as prl
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import os

# NOTE: check version
print(prl.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# SECTION: load THERMODB
# =======================================
# current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# thermodb reference
thermodb_reference = """
REFERENCES:
  CUSTOM-REF-1:
    TABLES:
      TABLE-ID: 103
      DESCRIPTION:
        "Thermodynamic properties of carbon dioxide (CO2)"
      DATA: []
      STRUCTURE:
        COLUMNS: [No., Name, Formula, State, Critical-Pressure, Critical-Temperature, Acentric-Factor, Enthalpy-of-Formation-Gas, Enthalpy-of-Formation-Liquid, Gibbs-Energy-of-Formation-Gas, Gibbs-Energy-of-Formation-Liquid]
        SYMBOL:  [None, None, None, None, Pc, Tc, AcFa, EnFo_IG, EnFo_LIQ, GiEnFo_IG, GiEnFo_LIQ]
        UNIT:    [None, None, None, None, MPa, K, None, kJ/mol, kJ/mol, kJ/mol, kJ/mol]
        CONVERSION: [None, None, None, None, 1, 1, 1, 1, 1, 1, 1]
      VALUES:
        - [1, "Carbon Dioxide", "CO2", "g", 7.38, 304.2, 0.225, -393.509, -393.509, -394.359, -394.359]
"""

# load thermodb
CO2 = None

# log
print(f'CO2: {CO2.check()}')


# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('CO2-g', CO2)


# NOTE: add thermodb rule
thermodb_config_file = os.path.join(current_dir, 'thermodb_config_link.yml')

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# SECTION: REACTION SYSTEM
# =======================================
# NOTE: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}
