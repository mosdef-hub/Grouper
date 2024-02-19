
<a name="readme-top"></a>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/kierannp/molGrouper">
    <img src="images/grouper.jpeg" alt="Logo" width="700" height="300">
  </a>

<h3 align="center">molGrouper</h3>

  <p align="center">
    A software package for creating and manipulating graphs of molecular groups.
    <br />
    <a href="https://github.com/kierannp/molGrouper"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/kierannp/molGrouper">View Demo</a>
    ·
    <a href="https://github.com/kierannp/molGrouper/issues">Report Bug</a>
    ·
    <a href="https://github.com/kierannp/molGrouper/issues">Request Feature</a>
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

* Need to install [`nauty`](https://pallini.di.uniroma1.it/) in packages in the base directory of `molGrouper`

### Installation

1. Clone the repo
```sh
git clone https://github.com/kierannp/molGrouper
cd molGrouper
```
2. Install with `pip`
```python
pip install .
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage

### Group graph initalization
```python
node_types = {
    'N': ['N1'], # amine
    'CO': ['C1', 'C2'], # carbonyl
    'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
    'C': ['C1'], # alkane
}


groupG = GroupGraph(node_types)
groupG.add_node('node1', 'N')
groupG.add_node('node2', 'CC')

groupG.add_edge('node1', 'N1', 'node2', 'C12')
```

### Exhaustive molecular space generation
```python
node_types = {
    'N': ['N1'], # amine
    'CO': ['C1', 'C2'], # carbonyl
    'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
    'C': ['C1', "C2", "C3", "C4"], # alkane
}

group_graph_space = generate_group_graph_space(
  n_nodes = 4, # number of groups to combine
  node_types = node_types
)
```

### Group graph to molecular graph
```python
node_type_to_smiles = {
    'N': 'N',
    'CO': 'C=O',
    'CC': 'C=C',
    'C': 'C',
}
node_port_to_atom_index = { # atom index is the atom that the port originates from
    'N': {'N1': 0}, 
    'CO': {'C1': 0, 'C2': 0},
    'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
    'C': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
}

mol_G = g.to_molecular_graph(
  node_type_to_smiles, 
  node_port_to_atom_index
)
```

### Molecular graph to SMILES
```python
from pysmiles import write_smiles

write_smiles(mol_G)
```


### Group graph from `mBuild.Compound`
```python
import mbuild as mb
from group_selfies import Group
from molGrouper import GroupGraph

mol = mb.load('CCCCCCCC', smiles=True) # octane molecule
groups = [Group('c3', 'C([H])([H])([H])(*1)'), Group('c2', 'C([H])([H])(*1)(*1)')]
groupG = GroupGraph()
groupG = groupG.from_mbuild(mol, groups)
```
### Conversion to `torch_geometric.Data`
```python
import molGrouper

node_types = {
    'CH2': ['C1', 'C2'], #methylene
    'CONH': ['C', 'N'], #amide
}

groupG = molGrouper.GroupGraph(node_types)
groupG.add_node('node1', 'CH2')
groupG.add_node('node2', 'CONH')
groupG.add_edge('node1', 'C1', 'node2', 'C')

group_featurizer = lambda node: torch.tensor([1, 0]) # dummy group featurizer

data = groupG.to_PyG_Data(group_featurizer, max_n_attachments=2)

data
# Data(x=[2,2], edge_index=[2,1], edge_attr=[1,4])
```



<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap
[X] Generate possible graphs with attachments points as limiter for number of connections 

[X] Process output of generation into python data structure

[X] Substiute nodes based on max attachments for graphs

[ ] Rdkit demo

[ ] group graph visualizaiton

[ ] Verify structures with synthesizability checks

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

Project Link: [https://github.com/kierannp/molGrouper](https://github.com/kierannp/molGrouper)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
[NSF GRFP](https://www.nsfgrfp.org/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/kierannp/molGrouper.svg?style=for-the-badge
[contributors-url]: https://github.com/kierannp/molGrouper/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/kierannp/molGrouper.svg?style=for-the-badge
[forks-url]: https://github.com/kierannp/molGrouper/network/members
[stars-shield]: https://img.shields.io/github/stars/kierannp/molGrouper.svg?style=for-the-badge
[stars-url]: https://github.com/kierannp/molGrouper/stargazers
[issues-shield]: https://img.shields.io/github/issues/kierannp/molGrouper.svg?style=for-the-badge
[issues-url]: https://github.com/kierannp/molGrouper/issues
[license-shield]: https://img.shields.io/badge/License-MIT-yellow.svg
[license-url]: https://github.com/kierannp/molGrouper/blob/master/LICENSE.txt