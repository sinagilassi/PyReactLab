# import packages/modules

# local
from .configs import __description__, __version__, packageName, \
    packageShortName


def intro():
    '''
    Package description
    '''
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def build_reaction_system(input_file, thermodb):
    '''
    Build a reaction system

    Parameters
    ----------
    input_file: str
        path of input file
    thermodb: dict
        thermodynamic database of components

    Returns
    -------
    ReactionSys
        reaction system
    '''
