# import libs
import pyreactlab as prl
from pyreactlab import ReactionCore
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pyreactlab.models import Reaction
from pythermodb_settings.models import Temperature, Pressure
import os
# locals
from tests.msource import (
    CO2,
    H2,
    CO,
    H2O,
    CH3OH,
    model_source
)

# NOTE: check version
print(prl.__version__)
print(ptdblink.__version__)

# =======================================
# SECTION: REACTION SYSTEM
# =======================================
# NOTE: reaction list format
rxn_1 = Reaction(
    name='Methanol Formation by CO2-Hydrogenation',
    reaction='CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)',
    components=[CO2, H2, CH3OH, H2O],
)

rxn_2 = Reaction(
    name='Reverse-Water-Gas-Shift',
    reaction='CO2(g) + H2(g) <=> CO(g) + H2O(g)',
    components=[CO2, H2, CO, H2O],
)

rxn_3 = Reaction(
    name='Methanol Formation by CO-Hydrogenation',
    reaction='CO(g) + 2H2(g) <=> CH3OH(g)',
    components=[CO, H2, CH3OH],
)

reactions = [rxn_1, rxn_2, rxn_3]

# =======================================
# SECTION: load THERMODB
# =======================================
# current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# =======================================
# SECTION: REACTION SYSTEM
# =======================================
# NOTE: summary
summary = prl.summary()
print(summary)

# NOTE: create reaction system
reaction_system = prl.create_gas_rxn(
    system_name='Methanol Synthesis',
    reactions=reactions,
    model_source=model_source,
    component_key="Formula-State"
)

# =======================================
# SECTION: REACTION SYSTEM PROPERTIES
# =======================================

# NOTE: Keq at T
res_ = reaction_system.reaction_equilibrium_constant(
    'Methanol Formation by CO2-Hydrogenation',
    temperature=Temperature(value=300.0, unit="K"),
    message="K_eq at 300 K",
)
print(f'K_eq: ')
print(res_)

# NOTE: select reaction
R1: ReactionCore = reaction_system.select_reaction(
    'Methanol Formation by CO2-Hydrogenation'
)
R2: ReactionCore = reaction_system.select_reaction(
    'Reverse-Water-Gas-Shift'
)
R3: ReactionCore = reaction_system.select_reaction(
    'Methanol Formation by CO-Hydrogenation'
)

# NOTE: equilibrium constant
# at std
print(f'K_eq: ')
print({R1.equilibrium_constant_std})

# at 300 K
# ! van't Hoff method
K_eq_300 = R1.cal_equilibrium_constant(
    Temperature(value=300.0, unit="K"), method="van't Hoff"
)
print(f'K_eq_300: ')
print(K_eq_300)

# ! shortcut van't Hoff method
K_eq_300 = R1.cal_equilibrium_constant(
    Temperature(value=300.0, unit="K"),
    method="shortcut van't Hoff"
)
print(f'K_eq_300: ')
print(K_eq_300)

# NOTE: energy of reaction
res_ = R1.cal_reaction_energy(
    Temperature(value=300.0, unit="K")
)
print(f'En_rxn: ')
print(res_)

# NOTE: formatting energies
res_ = reaction_system.component_formation_energies(
    component_name='CO2-g',
    temperature=Temperature(value=300.0, unit="K")
)
print(f'formation energies: ')
print(res_)

# =======================================
# SECTION: EQUILIBRIUM CALCULATION
# =======================================
# NOTE: mole fraction
mole_fraction = {
    'CO2-g': 1,
    'H2-g': 3,
    'CO-g': 1,
    'H2O-g': 0.1,
    'CH3OH-g': 0.3
}

mole = {
    'CO2-g': 1,
    'H2-g': 3,
    'CO-g': 1,
    'H2O-g': 0,
    'CH3OH-g': 0
}

# NOTE: activity inputs
# option 1: add activity inputs to datasource
# option 2: add activity inputs to inputs

# NOTE: input
inputs = {
    'mole': mole,
    'temperature': [100, "C"],
    'pressure': [1.0, "bar"],
}

# equilibrium calculation
res_ = reaction_system.equilibrium(
    inputs=inputs,
    conversion=['CO2-g'],
    method='least_squares',
    gas_mixture='ideal',
    solution='ideal',)
print(f'Equilibrium: {res_}')
