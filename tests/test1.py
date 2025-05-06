# import libs
import PyReactLab as prl
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
CO2 = ptdb.load_thermodb(f'{data_folder}/carbon dioxide.pkl')
H2 = ptdb.load_thermodb(f'{data_folder}/hydrogen.pkl')
CO = ptdb.load_thermodb(f'{data_folder}/carbon monoxide.pkl')
H2O = ptdb.load_thermodb(f'{data_folder}/water.pkl')
CH3OH = ptdb.load_thermodb(f'{data_folder}/methanol.pkl')
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
thub1.add_thermodb(CO2, name='CO2')
thub1.add_thermodb(H2, name='H2')
thub1.add_thermodb(CO, name='CO')
thub1.add_thermodb(H2O, name='H2O')
thub1.add_thermodb(CH3OH, name='CH3OH')

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

# NOTE: create reaction system
reaction_system = prl.create_rxn(
    system_name='Methanol Synthesis',
    reactions=reactions,
    model_source=model_source
)


# NOTE: mole fraction
mole_fraction = {
    'CO2': 0.1,
    'H2': 0.3,
    'CO': 0.2,
    'H2O': 0.1,
    'CH3OH': 0.3
}

# NOTE: model input
model_inputs = {
    'mole_fraction': mole_fraction,
    'temperature': [298.15, "K"],
    'pressure': [1.0, "bar"],
}
