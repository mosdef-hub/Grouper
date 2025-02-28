
<a name="readme-top"></a>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![CI](https://github.com/kierannp/Grouper/actions/workflows/CI.yaml/badge.svg)](https://github.com/kierannp/Grouper/actions/workflows/CI.yaml)



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/kierannp/Grouper">
    <img src="images/grouper.jpeg" alt="Logo" width="700" height="300">
  </a>

<h3 align="center">Grouper</h3>

  <p align="center">
    A software package for creating and manipulating graphs of molecular groups.
    <br />
    <a href="https://github.com/kierannp/Grouper"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/kierannp/Grouper">View Demo</a>
    ·
    <a href="https://github.com/kierannp/Grouper/issues">Report Bug</a>
    ·
    <a href="https://github.com/kierannp/Grouper/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

The fundamental data structure behind this package is based on a port graph, look [here](https://doi.org/10.1017/S0960129518000270) for an excellent description of the data structure.


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

*  Download and install [`nauty`](https://pallini.di.uniroma1.it/) according to their installation instructions in `Grouper/packages`

### Installation

1. Clone the repo
```sh
conda create -n Grouper python=3.9
conda activate Grouper
conda install librdkit-dev librdkit rdkit-dev rdkit-postgresql boost cmake rdkit eigen pybind11 openmp
conda install -c conda-forge nauty nlohmann_json

git clone https://github.com/kierannp/Grouper
cd Grouper
```
2. Install with `pip`
```python
python setup.py build_ext --inplace
python setup.py install
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage

### Group graph initalization
```python
from Grouper import GroupGraph

gG = GroupGraph()

# Adding nodes
gG.add_node(type = 'nitrogen', pattern = 'N', hubs = [0,0,0], is_smarts=False) # default of is_smarts is False
gG.add_node('nitrogen') # Once the type of the node has been specified we can use it again
gG.add_node(type = '', pattern = '[N]', hubs = [0,0,0], is_smarts=True) # Alternatively we can just use smarts

# Adding edges
gG.add_edge(src = (0,0), dst = (1,0), order=1) # In the format ((nodeID, srcPort), (nodeID, dstPort), bondOrder)
gG.add_edge(src = (1,1), dst = (2,0))
gG.add_edge(src = (2,1), dst = (0,1))

"""
Will make 
      N
     / \
    N - N
"""
```

### Group graph to SMILES
```python
smiles = gG.to_smiles()
```
### GroupGraph to AtomGraph
```python
atomG = gG.to_atomic_graph()
```



### Exhaustive chemical space generation
```python
import Grouper
from Grouper import Group, exhaustive_generate

node_defs = set()
# Define out node types that we will use to built our chemistries
node_defs.add(Group('nitrogen', 'N', [0,0,0]))
node_defs.add(Group('carbon', 'C', [0,0,0,0]))
node_defs.add(Group('oxygen', 'O', [0,0]))
node_defs.add(Group('benzene', 'c1ccccc1', [0,1,2,3,4,5]))

# Call method to enumerate possibilities
exhausted_space = exhaustive_generate(
    n_nodes = 4, 
    node_defs = node_defs, 
    input_file_path = '',
    nauty_path = '/path/to/nauty_X_X_X',
    num_procs = -1, # -1 utilizes all availible CPUs
)
```

### Conversion to `torch_geometric.Data`
```python
import Grouper
from Grouper.utils import convert_to_nx

def node_descriptor_generator(node_smiles):
    mol = rdkit.Chem.MolFromSmiles(node_smiles)
    desc = Descriptors.CalcMolDescriptors(mol)
    desc = {k: desc[k] for k in desc.keys() if (not isnan(desc[k]) or desc[k] is not None)}
    desc = [v for k,v in desc.items()] # flatten descriptors into single vector
    desc = torch.tensor(desc, dtype=torch.float64)
    return desc

nxG = convert_to_nx(gG)
max_ports = max(len(n.hubs) for n in node_defs)

data = nxG.to_PyG_Data(node_descriptor_generator, max_ports)
data.y = torch.tensor([rdkit.Chem.Descriptors.MolLogP(rdkit.Chem.MolFromSmiles(d))]) # here we utilize rdkit to estimate logP, but obviously can be generated another way 
```



<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the project
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Test modified project (`pytest`)
4. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
5. Push to the branch (`git push origin feature/AmazingFeature`)
6. Open a pull request

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Kieran Nehil-Puleo - nehilkieran@gmail.com

Cal Craven

Project Link: [https://github.com/kierannp/Grouper](https://github.com/kierannp/Grouper)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
[NSF GRFP](https://www.nsfgrfp.org/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/kierannp/Grouper.svg?style=for-the-badge
[contributors-url]: https://github.com/kierannp/Grouper/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/kierannp/Grouper.svg?style=for-the-badge
[forks-url]: https://github.com/kierannp/Grouper/network/members
[stars-shield]: https://img.shields.io/github/stars/kierannp/Grouper.svg?style=for-the-badge
[stars-url]: https://github.com/kierannp/Grouper/stargazers
[issues-shield]: https://img.shields.io/github/issues/kierannp/Grouper.svg?style=for-the-badge
[issues-url]: https://github.com/kierannp/Grouper/issues
[license-shield]: https://img.shields.io/badge/License-MIT-yellow.svg
[license-url]: https://github.com/kierannp/Grouper/blob/master/LICENSE.txt
