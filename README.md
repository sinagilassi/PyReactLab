# PyReactLab

![PyReactLab](https://drive.google.com/uc?export=view&id=1lJqEBwChS8iCN3WcfYVk6Sv-HqAlzEIB)

[![PyPI Downloads](https://static.pepy.tech/badge/mozichem-hub/month)](https://pepy.tech/projects/PyReactLab)
![PyPI](https://img.shields.io/pypi/v/PyReactLab)
![Python Version](https://img.shields.io/pypi/pyversions/PyReactLab.svg)
![License](https://img.shields.io/pypi/l/PyReactLab)
![Read the Docs](https://img.shields.io/readthedocs/PyReactLab)

PyReactLab is an open-source Python package designed to analyze chemical reactions and theoretically assess their thermodynamic feasibility. It provides a robust and flexible framework for studying chemical reaction systems, making it an ideal tool for researchers, engineers, and scientists working in chemistry, chemical engineering, and related disciplines.

## 🔎 Introduction

This package is used to analyze different chemical reactions and to thermodynamically check their feasibility and behavior under various conditions. PyReactLab is aimed at supporting theoretical studies, reaction pathway assessments, and reaction property evaluations without involving reactor design. It is especially useful in early-stage research and educational applications where understanding reaction viability is essential.

## ✨ Key Features

- **🧪 Reaction Feasibility Check**: Determine whether a given chemical reaction is thermodynamically feasible under specified conditions.

- **🔥 Thermodynamic Property Evaluation**: Calculate reaction enthalpies, Gibbs free energy changes, and equilibrium constants.

- **🔄 Multi-Reaction System Support**: Analyze sets of reactions, including coupled and parallel reactions.

- **📊 Flexible Input**: Support for user-defined reactions, species properties, and thermodynamic data.

- **🛠 Integration Ready**: Easily integrate PyReactLab into larger computational workflows or couple it with other modeling tools.

## 🚀 Installation

You can install PyReactLab using pip:

```bash
pip install PyReactLab
```

For the latest development version, you can install directly from GitHub:

```bash
pip install git+https://github.com/sinagilassi/PyReactLab.git
```

## 📚 Usage Example

Below is an example of how to use PyReactLab to analyze a reaction system:

```python
import PyReactLab as prl
from PyReactLab import Reaction
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import os

# Define reactions
reactions = [
    {'name': 'Methanol Formation by CO2-Hydrogenation', 'reaction': 'CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)'},
    {'name': 'Reverse-Water-Gas-Shift', 'reaction': 'CO2(g) + H2(g) <=> CO(g) + H2O(g)'},
    {'name': 'Methanol Formation by CO-Hydrogenation', 'reaction': 'CO(g) + 2H2(g) <=> CH3OH(g)'}
]

# Load thermodynamic data
current_dir = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(current_dir, 'data')
CO2 = ptdb.load_thermodb(f'{data_folder}/carbon dioxide-1.pkl')
H2 = ptdb.load_thermodb(f'{data_folder}/hydrogen-1.pkl')
CO = ptdb.load_thermodb(f'{data_folder}/carbon monoxide-1.pkl')
H2O = ptdb.load_thermodb(f'{data_folder}/water-1.pkl')
CH3OH = ptdb.load_thermodb(f'{data_folder}/methanol-1.pkl')

# Configure thermodynamic database
thub1 = ptdblink.init()
thub1.add_thermodb('CO2-g', CO2)
thub1.add_thermodb('H2-g', H2)
thub1.add_thermodb('CO-g', CO)
thub1.add_thermodb('H2O-g', H2O)
thub1.add_thermodb('CH3OH-g', CH3OH)
thermodb_config_file = os.path.join(current_dir, 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)
datasource, equationsource = thub1.build()

# Create reaction system
model_source = {"datasource": datasource, "equationsource": equationsource}

# create a reaction system
# option 1: create_gas_rxn
# option 2: create_liquid_rxn
# option 3: create_rxn
reaction_system = prl.create_rxn(
    system_name='Methanol Synthesis',
    reactions=reactions,
    model_source=model_source
)

# Calculate equilibrium constant at 300 K
res_ = reaction_system.reaction_equilibrium_constant(
    'Methanol Formation by CO2-Hydrogenation',
    [300.0, "K"],
    message="K_eq at 300 K",
)
print(f'K_eq: {res_}')

# Perform equilibrium calculation
inputs = {
    'mole': {'CO2-g': 1, 'H2-g': 3, 'CO-g': 1, 'H2O-g': 0.001, 'CH3OH-g': 0.001},
    'temperature': [100, "C"],
    'pressure': [1.0, "bar"],
}
res_ = reaction_system.equilibrium(
    inputs=inputs,
    conversion=['CO2-g'],
    method='minimize',
    gas_mixture='ideal',
    solution='ideal',
)
print(f'Equilibrium: {res_}')
```

This example demonstrates how to define reactions, load thermodynamic data, configure the database, and perform equilibrium calculations.

The result are as:

```python
{'feed':
    {'mole': {
        'value': {
            'H2': 3,
            'CH3OH': 0.001,
            'CO': 1,
            'CO2': 1,
            'H2O': 0.001
            },
            'unit': 'mol'
            },
    'mole_total': {
        'value': 5.002,
        'unit': 'mol'
        },
    'mole_fraction': {
        'value':{
            'H2': 0.5997600959616153,
            'CH3OH': 0.0001999200319872051,
            'CO': 0.1999200319872051,
            'CO2': 0.1999200319872051,
            'H2O': 0.0001999200319872051
            },
            'unit': 'dimensionless'
            },
    'mole_fraction_sum': {
        'value': 0.9999999999999999,
        'unit': 'dimensionless'
        }
    },
'equilibrium':
    {'mole':
        {'value':{
            'H2': 1.5659480019741632,
            'CH3OH': 0.7177747496801842,
            'CO': 0.2837277489852842,
            'CO2': 0.9994975013345316,
            'H2O': 0.001502498665468366
            },
            'unit': 'mol'
            },
        'mole_total': {
            'value': 3.5684505006396314,
            'unit': 'mol'
            },
        'mole_fraction': {
            'value':{
                'H2': 0.4388313643956874,
                'CH3OH': 0.20114465635757753,
                'CO': 0.07951006996858358,
                'CO2': 0.2800928585545394,
                'H2O': 0.0004210507236121249
                },
                'unit': 'dimensionless'
                },
        'mole_fraction_sum': {
            'value': 0.9999999999999999,
            'unit': 'dimensionless'
            }
    },
'extent_of_reaction': {
    'value': array([2.48905604e-05, 4.77608105e-04, 7.16749859e-01]),
    'unit':'dimensionless'
},
'temperature': {
    'value': 373.15,
    'unit': 'K'
    },
'pressure': {
    'value': 1.0,
    'unit': 'bar'
    },
'optimization_fun': {
    'value': 3.2676427381965915e-11,
    'unit':'dimensionless'
    },
'conversion': {
    'CO2': {
        'value': 0.0502498665468365,
        'unit': '%'
        }
        },
'gas_mixture': 'ideal',
'solution': 'ideal',
'computation_time': {
    'value': 1.5210797786712646,
    'unit': 's'}
    }
```

## 🤝 Contributions

Contributions, issues, and feature requests are welcome!
Feel free to check the issues page or submit a pull request.

## 📝 License

This project is licensed under the MIT License. You are free to use, modify, and distribute this software in your own applications or projects. However, if you choose to use this app in another app or software, please ensure that my name, Sina Gilassi, remains credited as the original author. This includes retaining any references to the original repository or documentation where applicable. By doing so, you help acknowledge the effort and time invested in creating this project.

## ❓ FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/)
