# import libs
import re
from math import exp
from typing import Dict, Any, List, Optional, Literal
from ..configs import (
    R_CONST_J__molK, DATASOURCE, EQUATIONSOURCE,
    PRESSURE_REF_Pa, TEMPERATURE_REF_K
)


class ChemReactUtils:
    """A utility class for chemical reaction analysis."""

    # variables
    _system_inputs = None
    # universal gas constant [J/mol.K]
    R = R_CONST_J__molK
    # temperature [K]
    T_Ref = TEMPERATURE_REF_K
    # pressure [bar]
    P_Ref = PRESSURE_REF_Pa/1e5

    def __init__(self,
                 reaction_mode_symbol: Literal["<=>"] = "<=>"):
        """
        Initialize the ChemReactUtils class.

        Parameters
        ----------
        reaction_mode_symbol : str, optional
            The symbol used to separate reactants and products in a reaction equation.
        """
        self.reaction_mode_symbol = reaction_mode_symbol

    @property
    def system_inputs(self) -> Dict[str, Any]:
        """Get the system inputs."""
        return self._system_inputs

    def count_carbon(self, molecule: str, coefficient: float) -> float:
        """
        Count the number of carbon atoms in a molecule.

        Parameters
        ----------
        molecule : str
            The chemical formula of the molecule.
        coefficient : float
            The coefficient of the molecule in the reaction.

        Returns
        -------
        float
            The number of carbon atoms in the molecule multiplied by the coefficient.
        """
        try:
            # NOTE: check molecule
            if not isinstance(molecule, str):
                raise ValueError("Molecule must be a string.")

            # NOTE: check coefficient
            if not isinstance(coefficient, (int, float)):
                raise ValueError("Coefficient must be an integer or float.")

            # NOTE: Check if the molecule contains carbon atoms
            if re.search(r'C(?![a-z])', molecule):
                carbon_count = len(re.findall(
                    r'C(?![a-z])', molecule)) * coefficient
                return carbon_count
            else:
                return 0.0
        except Exception as e:
            raise Exception(
                f"Error counting carbon in molecule '{molecule}': {e}")

    def analyze_reaction(
            self,
            reaction_pack: Dict[str, str]) -> Dict[str, Any]:
        """
        Analyze a chemical reaction and extract relevant information.

        Parameters
        ----------
        reaction_pack : dict
            A dictionary containing the reaction and its name.

        Returns
        -------
        dict
            A dictionary containing the analyzed reaction data, including reactants,
            products, reaction coefficient, and carbon count.
        """
        try:
            # NOTE: check reaction_pack
            if not isinstance(reaction_pack, dict):
                raise ValueError("reaction_pack must be a dictionary.")

            if 'reaction' not in reaction_pack or 'name' not in reaction_pack:
                raise ValueError(
                    "reaction_pack must contain 'reaction' and 'name' keys.")

            # SECTION: extract data from reaction
            reaction = reaction_pack['reaction']
            name = reaction_pack['name']

            # Split the reaction into left and right sides
            sides = reaction.split(self.reaction_mode_symbol.strip())

            # Define a regex pattern to match reactants/products
            # pattern = r'(\d*)?(\w+)\((\w)\)'
            pattern = r'(\d*\.?\d+)?(\w+)\((\w)\)'

            # Extract reactants
            reactants = re.findall(pattern, sides[0])
            reactants = [{'coefficient': float(r[0]) if r[0] else float(
                1), 'molecule': r[1], 'state': r[2]} for r in reactants]

            # Extract products
            products = re.findall(pattern, sides[1])
            products = [{'coefficient': float(p[0]) if p[0] else float(
                1), 'molecule': p[1], 'state': p[2]} for p in products]

            # reaction coefficient
            reaction_coefficient = 0
            for item in reactants:
                reaction_coefficient += item['coefficient']
            for item in products:
                reaction_coefficient -= item['coefficient']

            # SECTION: Carbon count for each component
            carbon_count = {}
            for r in reactants:
                carbon_count[r['molecule']] = self.count_carbon(
                    r['molecule'], r['coefficient'])
            for p in products:
                carbon_count[p['molecule']] = self.count_carbon(
                    p['molecule'], p['coefficient'])

            # res
            res = {
                'name': name,
                'reaction': reaction,
                'reactants': reactants,
                'products': products,
                'reaction_coefficient': reaction_coefficient,
                'carbon_count': carbon_count
            }

            return res
        except Exception as e:
            raise Exception(f"Error analyzing reaction '{reaction}': {e}")

    def analyze_overall_reactions(
            self,
            reactions: List[Dict[str, str]]) -> Dict[str, List[str]]:
        """
        Analyze a list of chemical reactions and classify species as consumed, produced, or intermediate.

        Parameters
        ----------
        reactions : list
            A list of dictionaries, each containing a reaction string and its name.

        Returns
        -------
        dict
            A dictionary containing three lists: 'consumed', 'produced', and 'intermediate'.
            - 'consumed': List of species consumed in the reactions.
            - 'produced': List of species produced in the reactions.
            - 'intermediate': List of species that are both consumed and produced in the reactions.
        """
        try:
            # Initialize sets for all reactants and products
            all_reactants = set()
            all_products = set()

            # Iterate over reactions
            for reaction in reactions:
                # Split the reaction into left and right sides
                sides = reaction['reaction'].split(
                    self.reaction_mode_symbol.strip())

                # Define a regex pattern to match reactants/products
                pattern = r'(\d*)?(\w+)\((\w)\)'

                # Extract reactants
                reactants = re.findall(pattern, sides[0])
                reactants = [r[1] for r in reactants]

                # Extract products
                products = re.findall(pattern, sides[1])
                products = [p[1] for p in products]

                # Update sets
                all_reactants.update(reactants)
                all_products.update(products)

            # Classify species
            consumed = list(all_reactants - all_products)
            produced = list(all_products - all_reactants)
            intermediate = list(all_reactants & all_products)

            # res
            res = {
                'consumed': consumed,
                'produced': produced,
                'intermediate': intermediate
            }

            return res
        except Exception as e:
            raise Exception(f"Error analyzing overall reactions: {e}")

    def define_component_id(self, reaction_res):
        '''
        Define component ID

        Parameters
        ----------
        reaction_res: dict
            reaction_res

        Returns
        -------
        component_list: list
            component list
        component_dict: dict
            component dict
        comp_list: list
            component list
        comp_coeff: list
            component coefficient
        '''
        try:
            # NOTE: component list
            component_list = []

            # Iterate over reactions and extract reactants and products
            for item in reaction_res:
                for reactant in reaction_res[item]['reactants']:
                    component_list.append(reactant['molecule'])
                for product in reaction_res[item]['products']:
                    component_list.append(product['molecule'])

            # remove duplicate
            component_list = list(set(component_list))

            # component id: key, value
            component_dict = {}
            for i, item in enumerate(component_list):
                component_dict[item] = i

            # Initialize the component list
            comp_list = [{i: 0.0 for i in component_dict.keys()}
                         for _ in range(len(reaction_res))]

            # NOTE: Iterate over reactions and components
            for j, reaction in enumerate(reaction_res):
                for item in component_dict.keys():
                    # Check reactants
                    for reactant in reaction_res[reaction]['reactants']:
                        if reactant['molecule'] == item:
                            comp_list[j][item] = -1 * \
                                float(reactant['coefficient'])

                    # Check products
                    for product in reaction_res[reaction]['products']:
                        if product['molecule'] == item:
                            comp_list[j][item] = float(product['coefficient'])

            # Convert comp_list to comp_matrix
            comp_coeff = [[comp_list[j][item] for item in component_dict.keys()]
                          for j in range(len(reaction_res))]

            # res
            return component_list, component_dict, comp_list, comp_coeff
        except Exception as e:
            raise Exception(f"Error defining component ID: {e}")
