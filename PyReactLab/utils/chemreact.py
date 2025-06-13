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

    # # NOTE: variables
    # system inputs
    _system_inputs = None
    # universal gas constant [J/mol.K]
    R = R_CONST_J__molK
    # temperature [K]
    T_Ref = TEMPERATURE_REF_K
    # pressure [bar]
    P_Ref = PRESSURE_REF_Pa/1e5

    # available phases
    available_phases = ['gas', 'liquid', 'aqueous', 'solid']

    def __init__(
        self,
        reaction_mode_symbol: Literal["<=>"] = "<=>"
    ):
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
        # check
        if self._system_inputs is None:
            raise ValueError("System inputs are not set.")
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
            reaction_pack: Dict[str, str],
            phase_rule: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Analyze a chemical reaction and extract relevant information.

        Parameters
        ----------
        reaction_pack : dict
            A dictionary containing the reaction and its name.
        phase_rule : str, optional
            The phase of the reaction, which can be 'gas', 'liquid', 'aqueous', or 'solid'.

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

            # NOTE: check phase
            # set phase
            phase_set = self.phase_rule_analysis(phase_rule)

            # SECTION: extract data from reaction
            reaction = reaction_pack['reaction']
            name = reaction_pack['name']

            # Split the reaction into left and right sides
            sides = reaction.split(self.reaction_mode_symbol.strip())

            # Define a regex pattern to match reactants/products
            # pattern = r'(\d*)?(\w+)\((\w)\)'
            # pattern = r'(\d*\.?\d+)?(\w+)\((\w)\)'
            pattern = r'(?:(\d*\.?\d+)\s*)?([A-Z][a-zA-Z0-9]*)\s*(?:\((\w)\))?'

            # Extract reactants
            reactants = re.findall(pattern, sides[0])
            reactants = [
                {
                    'coefficient': float(r[0]) if r[0] else float(1),
                    'molecule': r[1],
                    'state': r[2] if r[2] else phase_set
                } for r in reactants
            ]

            # NOTE: reactants full name
            reactants_names = []
            # loop over reactants
            for i, item in enumerate(reactants):
                # ! check phase_set and phase_rule
                if phase_rule is None:
                    # check item state
                    if item['state'] == 'empty':
                        raise ValueError(
                            f"Phase rule is empty but reactant '{item['molecule']}' has state '{item['state']}'.")
                else:
                    # check item state
                    if item['state'] != phase_set:
                        raise ValueError(
                            f"Phase rule is '{phase_set}' but reactant '{item['molecule']}' has state '{item['state']}'.")

                # generate full name
                full_name = item['molecule'] + "-" + item['state']
                # append to list
                reactants_names.append(full_name)
                # update source
                reactants[i]['molecule_state'] = full_name

            # Extract products
            products = re.findall(pattern, sides[1])
            products = [
                {
                    'coefficient': float(p[0]) if p[0] else float(1),
                    'molecule': p[1],
                    'state': p[2] if p[2] else phase_set
                } for p in products
            ]

            # NOTE: products full name
            products_names = []
            # loop over products
            for i, item in enumerate(products):
                # ! check phase_set and phase_rule
                if phase_rule is None:
                    # check item state
                    if item['state'] == 'empty':
                        raise ValueError(
                            f"Phase rule is empty but product '{item['molecule']}' has state '{item['state']}'.")
                else:
                    # check item state
                    if item['state'] != phase_set:
                        raise ValueError(
                            f"Phase rule is '{phase_set}' but product '{item['molecule']}' has state '{item['state']}'.")

                # generate full name
                full_name = item['molecule'] + "-" + item['state']
                # append to list
                products_names.append(full_name)
                # update source
                products[i]['molecule_state'] = full_name

            # reaction coefficient
            reaction_coefficient = 0
            for item in reactants:
                reaction_coefficient += item['coefficient']
            for item in products:
                reaction_coefficient -= item['coefficient']

            # SECTION: Carbon count for each component
            carbon_count = {}
            for r in reactants:
                carbon_count[r['molecule_state']] = self.count_carbon(
                    r['molecule'], r['coefficient'])
            for p in products:
                carbon_count[p['molecule_state']] = self.count_carbon(
                    p['molecule'], p['coefficient'])

            # SECTION: reaction state
            reaction_state = {}
            for r in reactants:
                # set
                reaction_state[r['molecule_state']] = r['state']
            for p in products:
                # set
                reaction_state[p['molecule_state']] = p['state']

            # NOTE: reaction phase
            # reaction
            reaction_phase = self.determine_reaction_phase(
                reaction_state)

            # NOTE: unique states
            state_count = self.count_reaction_states(
                reaction_state)

            # res
            res = {
                'name': name,
                'reaction': reaction,
                'reactants': reactants,
                'reactants_names': reactants_names,
                'products': products,
                'products_names': products_names,
                'reaction_coefficient': reaction_coefficient,
                'carbon_count': carbon_count,
                'reaction_state': reaction_state,
                'reaction_phase': reaction_phase,
                'state_count': state_count,
            }

            return res
        except Exception as e:
            raise Exception(f"Error analyzing reaction '{reaction}': {e}")

    def phase_rule_analysis(self, phase_rule: Optional[str] = None) -> str:
        """
        Analyze the phase rule of a reaction.

        Parameters
        ----------
        phase_rule : str, optional
            The phase rule of the reaction.

        Returns
        -------
        phase_symbol : str
            The phase symbol of the reaction, which can be 'g', 'l', or 'empty'.
        """
        try:
            # SECTION: check phase rule
            # Check if phase_rule is None
            if phase_rule is None or phase_rule == 'None':
                # set default phase
                return 'empty'

            # SECTION: check phase rule
            # Check if the phase rule is valid
            if phase_rule not in self.available_phases:
                raise ValueError(
                    f"Phase rule must be {', '.join(self.available_phases)}.")

            # check phase
            if phase_rule == 'gas':
                phase_symbol = 'g'
            elif phase_rule == 'liquid':
                phase_symbol = 'l'
            elif phase_rule == 'aqueous':
                phase_symbol = 'aq'
            elif phase_rule == 'solid':
                phase_symbol = 's'
            else:
                phase_symbol = 'empty'

            return phase_symbol
        except Exception as e:
            raise Exception(f"Error analyzing phase rule: {e}")

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
                # NOTE: Split the reaction into left and right sides
                sides = reaction['reaction'].split(
                    self.reaction_mode_symbol.strip())

                # NOTE: Define a regex pattern to match reactants/products
                # pattern = r'(\d*)?(\w+)\((\w)\)'
                pattern = r'(?:(\d*\.?\d+)\s*)?([A-Z][a-zA-Z0-9]*)\s*(?:\((\w)\))?'

                # SECTION: Extract reactants
                reactants = re.findall(pattern, sides[0])
                reactants = [r[1] for r in reactants]

                # SECTION: Extract products
                products = re.findall(pattern, sides[1])
                products = [p[1] for p in products]

                # NOTE: Update sets
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

    def analyze_overall_reactions_v2(
        self,
        reactions: Dict[str, Any]
    ) -> Dict[str, List[str]]:
        """
        Analyze a list of chemical reactions and classify species as consumed, produced, or intermediate (version 2).

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
            for reaction_name, reaction_value in reactions.items():
                # SECTION: Extract reactants
                reactants = reaction_value['reactants_names']

                # SECTION: Extract products
                products = reaction_value['products_names']

                # NOTE: Update sets
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

            # SECTION: Iterate over reactions and extract reactants and products
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

            # SECTION: Initialize the component list
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

    def define_component_id_v2(self, reaction_res):
        '''
        Define component ID (version 2)

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

            # SECTION: Iterate over reactions and extract reactants and products
            for item in reaction_res:
                # reactants
                for reactant in reaction_res[item]['reactants']:
                    component_list.append(reactant['molecule_state'])
                # products
                for product in reaction_res[item]['products']:
                    component_list.append(product['molecule_state'])

            # remove duplicate
            component_list = list(set(component_list))

            # component id: key, value
            component_dict = {}

            # loop over component list
            for i, item in enumerate(component_list):
                component_dict[item] = i

            # SECTION: Initialize the component list
            comp_list = [{i: 0.0 for i in component_dict.keys()}
                         for _ in range(len(reaction_res))]

            # NOTE: Iterate over reactions and components
            for j, reaction in enumerate(reaction_res):
                for item in component_dict.keys():
                    # Check reactants
                    for reactant in reaction_res[reaction]['reactants']:
                        if reactant['molecule_state'] == item:
                            comp_list[j][item] = -1 * \
                                float(reactant['coefficient'])

                    # Check products
                    for product in reaction_res[reaction]['products']:
                        if product['molecule_state'] == item:
                            comp_list[j][item] = float(product['coefficient'])

            # Convert comp_list to comp_matrix
            comp_coeff = [[comp_list[j][item] for item in component_dict.keys()]
                          for j in range(len(reaction_res))]

            # res
            return component_list, component_dict, comp_list, comp_coeff
        except Exception as e:
            raise Exception(f"Error defining component ID: {e}")

    def state_name_set(self, state_set: set) -> List[str]:
        '''
        Convert state set to full names

        Parameters
        ----------
        state_set: set
            Set of states

        Returns
        -------
        state_names: list
            List of full state names
        '''
        try:
            state_dict = {
                'g': 'gas',
                'l': 'liquid',
                'aq': 'aqueous',
                's': 'solid'
            }

            # Map the states from the set to their full names
            return [state_dict[state] for state in state_set]
        except Exception as e:
            raise Exception(f"Error converting state set to full names: {e}")

    def determine_reaction_phase(self, reaction_dict: Dict[str, str]) -> str:
        '''
        Determine the phase of a reaction based on the states of its components.

        Parameters
        ----------
        reaction_dict: dict
            A dictionary where keys are component names and values are their states.

        Returns
        -------
        str
            The phase of the reaction, which can be 'gas', 'liquid', 'aqueous', 'solid', or a combination of these.
        '''
        try:
            # Collect the states from the values in the dictionary
            available_states = set(reaction_dict.values())

            # Convert the states to full names
            state_names = self.state_name_set(available_states)

            # Determine phase based on the number of unique states
            if len(state_names) == 1:
                return f'{state_names[0]}'
            else:
                return f'{"-".join(state_names)}'
        except Exception as e:
            raise Exception(f"Error determining reaction phase: {e}")

    def count_reaction_states(self, reaction_dict: Dict[str, str]) -> Dict[str, int]:
        '''
        Counts the number of unique states in a reaction as g, l, aq, or s.

        Parameters
        ----------
        reaction_dict: dict
            A dictionary where keys are component names and values are their states.

        Returns
        -------
        dict
            A dictionary with the counts of each state (g, l, aq, s).
        '''
        try:
            # Collect the states from the values in the dictionary
            available_states = reaction_dict.values()

            # how many g, l, aq, or s
            state_count = {
                'g': 0,
                'l': 0,
                'aq': 0,
                's': 0
            }

            # Count the occurrences of each state
            for state in available_states:
                if state in state_count:
                    state_count[state] += 1

            return state_count

        except Exception as e:
            raise Exception(f"Error determining reaction phase: {e}")

    def reaction_phase_analysis(
        self,
        reaction_res: Dict[str, Any],
    ):
        '''
        Analyze the reaction phase and separate reactants and products by their phases.

        Parameters
        ----------
        reaction_res: dict
            A dictionary containing the reaction results, including reactants and products.

        Returns
        -------

        '''
        try:
            # NOTE: initialize phase dict
            phase_dict = {
                'g': [],
                'l': [],
                'aq': [],
                's': []
            }

            # SECTION: Iterate over reactions and classify reactants and products by phase
            for reaction_name, reaction_data in reaction_res.items():
                # reactants
                for reactant in reaction_data['reactants']:
                    phase = reactant['state']
                    if phase in phase_dict:
                        # ! molecule state
                        phase_dict[phase].append(reactant['molecule_state'])
                    else:
                        raise ValueError(
                            f"Unknown phase '{phase}' for reactant '{reactant['molecule']}'.")

                # products
                for product in reaction_data['products']:
                    phase = product['state']
                    if phase in phase_dict:
                        # ! molecule state
                        phase_dict[phase].append(product['molecule_state'])
                    else:
                        raise ValueError(
                            f"Unknown phase '{phase}' for product '{product['molecule']}'.")

            # NOTE: remove duplicates in each phase
            for phase in phase_dict:
                phase_dict[phase] = list(set(phase_dict[phase]))

            # res
            return phase_dict
        except Exception as e:
            raise Exception(f"Error analyzing reaction phase: {e}")
