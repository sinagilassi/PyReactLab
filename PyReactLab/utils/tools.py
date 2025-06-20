# TOOLS
# import libs
from ..configs import DATASOURCE, EQUATIONSOURCE


def model_source_checker(model_source: dict) -> bool:
    """
    Check if the model source is valid.

    Parameters
    ----------
    model_source : dict
        Model source parameters as a dictionary defined as:
        - `datasource`: Dictionary containing the model data.
        - `equationsource`: Dictionary containing the equations.

    Returns
    -------
    bool
        True if the model source is valid, False otherwise.
    """
    try:
        # Check if the model source is a dictionary and contains the required keys
        if (not isinstance(model_source, dict) or
                not all(key in model_source for key in [DATASOURCE, EQUATIONSOURCE])):
            raise ValueError(
                f"Model source must be a dictionary with {DATASOURCE} and {EQUATIONSOURCE} keys.")

        # Check if the datasource and equationsource are dictionaries
        if (not isinstance(model_source[DATASOURCE], dict) or
                not isinstance(model_source[EQUATIONSOURCE], dict)):
            raise ValueError(
                f"'{DATASOURCE}' and '{EQUATIONSOURCE}' must be dictionaries.")

        return True
    except Exception as e:
        raise RuntimeError(f"Model source validation failed: {e}") from e
