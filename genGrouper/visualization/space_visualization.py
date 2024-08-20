import typing as t
import matplotlib.pyplot as plt
from numpy.typing import ArrayLike
from grakel.kernels import Kernel



def plot_uncharacterized_space(space_set : t.Set[str], kernel: Kernel, cutoff: float, options: dict) -> plt.Figure:
    pass

def plot_characterized_space(space_set : t.Set[str], characterized_values: t.Dict[str, t.Union[float, ArrayLike, str]], kernel: Kernel, cutoff: float, options: dict) -> plt.Figure:
    pass

def plot_partiallly_characterized_space(space_set : t.Set[str], characterized_values: t.Dict[str, t.Union[float, ArrayLike, str]], kernel: Kernel, cutoff: float, options: dict) -> plt.Figure:
    pass

def plot_feature_similarity(feature_set: t.Mapping[str, ArrayLike], kernel: t.Callable[[ArrayLike, ArrayLike], float], cutoff: float, options: dict) -> plt.Figure:
    pass
