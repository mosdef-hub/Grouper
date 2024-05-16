
import typing as t

def multi_to_pair(multi: int, max_multi: int) -> t.Tuple[int, int]:
    """
    Convert a linear index `multi` to 2D coordinates (x, y) within a square grid.

    Parameters:
    - multi (int): Linear index to be converted.
    - max_multi (int): Maximum allowed value for the linear index.

    Returns:
    - Tuple[int, int]: Two-dimensional coordinates (x, y) within the grid.
    """
    if 1 <= multi <= max_multi:
        # Calculate x and y values based on the input multi
        x = (multi - 1) % int(max_multi**.5) + 1
        y = (multi - 1) // int(max_multi**.5) + 1

        x, y = x-1, y-1 # convert to 0-indexed
        return x, y
    else:
        raise ValueError("Input multi must be in the range 1 to max_multi.")