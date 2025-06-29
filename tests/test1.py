# import libs
import PyReactLab as prl
from PyReactLab import Reaction
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import os

# NOTE: check version
print(prl.__version__)
print(ptdblink.__version__)

# =======================================
# SECTION: REACTION SYSTEM
# =======================================
# NOTE: reaction list format
reaction_1 = 'CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)'
reaction_2 = 'CO2(g) + H2(g) <=> CO(g) + H2O(g)'
reaction_3 = 'CO(g) + 2H2(g) <=> CH3OH(g)'

reactions = [
    {
        'name': 'Methanol Formation by CO2-Hydrogenation',
        'reaction': reaction_1
    },
    {
        'name': 'Reverse-Water-Gas-Shift',
        'reaction': reaction_2
    },
    {
        'name': 'Methanol Formation by CO-Hydrogenation',
        'reaction': reaction_3
    }
]

# =======================================
# SECTION: load THERMODB
# =======================================
# current directory
current_dir = os.path.dirname(os.path.abspath(__file__))
# data folder
data_folder = os.path.join(current_dir, 'data')

# load thermodb
CO2 = ptdb.load_thermodb(f'{data_folder}/carbon dioxide-1.pkl')
H2 = ptdb.load_thermodb(f'{data_folder}/hydrogen-1.pkl')
CO = ptdb.load_thermodb(f'{data_folder}/carbon monoxide-1.pkl')
H2O = ptdb.load_thermodb(f'{data_folder}/water-1.pkl')
CH3OH = ptdb.load_thermodb(f'{data_folder}/methanol-1.pkl')
# log
print(f'CO2: {CO2.check()}')
print(f'H2: {H2.check()}')
print(f'CO: {CO.check()}')
print(f'H2O: {H2O.check()}')
print(f'CH3OH: {CH3OH.check()}')

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('CO2-g', CO2)
thub1.add_thermodb('H2-g', H2)
thub1.add_thermodb('CO-g', CO)
thub1.add_thermodb('H2O-g', H2O)
thub1.add_thermodb('CH3OH-g', CH3OH)

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

# NOTE: summary
summary = prl.summary()
print(summary)

# NOTE: create reaction system
reaction_system = prl.create_gas_rxn(
    system_name='Methanol Synthesis',
    reactions=reactions,
    model_source=model_source
)

# SECTION: REACTION SYSTEM PROPERTIES

# Keq at T
res_ = reaction_system.reaction_equilibrium_constant(
    'Methanol Formation by CO2-Hydrogenation',
    [300.0, "K"],
    message="K_eq at 300 K",
)
print(f'K_eq: {res_}')
# NOTE: select reaction
R1: Reaction = reaction_system.select_reaction(
    'Methanol Formation by CO2-Hydrogenation')
R2: Reaction = reaction_system.select_reaction('Reverse-Water-Gas-Shift')
R3: Reaction = reaction_system.select_reaction(
    'Methanol Formation by CO-Hydrogenation')

# NOTE: equilibrium constant
# at std
print(f'K_eq: {R1.equilibrium_constant_std}')
# at 300 K
K_eq_300 = R1.cal_equilibrium_constant([300.0, "K"])
print(f'K_eq_300: {K_eq_300}')
K_eq_300 = R1.cal_equilibrium_constant(
    [300.0, "K"], method="shortcut van't Hoff")
print(f'K_eq_300: {K_eq_300}')

# NOTE: energy of reaction
res_ = R1.cal_reaction_energy([300.0, "K"])
print(f'En_rxn: {res_}')

# NOTE: formatting energies
res_ = reaction_system.component_formation_energies(
    component_name='CO2-g',
    temperature=[300.0, "K"],
)
print(f'formation energies: {res_}')

# SECTION: EQUILIBRIUM CALCULATION

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
