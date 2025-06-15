from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyReactLab'
VERSION = '0.3.0'
AUTHOR = 'Sina Gilassi'
EMAIL = "<sina.gilassi@gmail.com>"
DESCRIPTION = 'PyReactLab is an open-source Python package designed to analyze chemical reactions and theoretically assess their thermodynamic feasibility.'
LONG_DESCRIPTION = "PyReactLab is an open-source Python package designed to perform comprehensive analysis of chemical reactions, with a primary focus on the theoretical assessment of their thermodynamic feasibility. It provides a versatile and extensible computational framework that allows users to model, analyze, and interpret chemical reaction systems under various conditions. The package is particularly useful for researchers and engineers in the field of chemical engineering, reaction kinetics, and thermodynamics. PyReactLab aims to facilitate the understanding of reaction mechanisms, equilibrium constants, and kinetic parameters, making it a valuable tool for both academic and industrial applications."

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    include_package_data=True,  # Make sure to include non-Python files
    # Add both config and data files
    package_data={'': ['config/*.yml', 'config/*.csv', 'data/*.csv',
                       'data/*.yml', 'plugin/*.yml', 'plugin/*.csv']},
    license='MIT',
    license_files=[],
    install_requires=['numpy',
                      'PyYAML',
                      'PyCUC',
                      'scipy',
                      'PyThermoModels',
                      'PyThermoDB'],
    extras_require={
        "plotting": ["matplotlib"],
    },
    keywords=[
        'python', 'chemical-engineering', 'process-modeling',
        'reaction-analysis', 'equilibrium-constant', 'kinetics',
        'thermodynamic-analysis', 'reaction-kinetics', 'thermodynamics',
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.13",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.10',
)
