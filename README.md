# netfibGen
Library to generate fibre networks used in computational homogenisation analysis. 

This library has been developed by me, Felipe Figueredo Rocha, during my PhD in Computational Modelling 2015 - 2019 at LNCC, Brazil. It has been adapted to be made public available in 2021. 

Principle: The network of fibres is computationally generated by providing a set of target properties. Initially, for a given average fibre orientation, say (measured with respect to the horizontal axis) and for a certain number of fibres, a homogeneous network is generated containing two families of fibres symmetrically oriented. Crossing-points are considered to be junctions, which are the extremes of computational segments composing the fibres. Then, the position of each node is individually perturbed in a random manner in terms of distance and direction, limited by a circle of radius. It turns out that the orientation of each fibre results from a combination of a given mean value and the perturbation. Once the network has been built, spatial distribution of material properties and fibre areas are selected. In the study cases presented below we consider either properties constant for all fibres, properties randomly sampled from a known probability density function (e.g. a normal distribution), or specifically modified in a specific region of the RVE, such as a bands or balls. For the sake of simplicity, in the forthcoming examples the  damage threshold stress and the fibre area are considered sources of heterogeneity.

Please cite at least one of the following documents if this library has been useful to you:

@article{Rocha2021, 
title = {Damage-driven strain localisation in networks of fibres: A computational homogenisation approach}, 
author = {Felipe Figueredo Rocha and Pablo Javier Blanco and Pablo Javier Sánchez and Eduardo {de Souza Neto} and Raúl Antonino Feijóo}, 
journal = {Computers & Structures}, 
volume = {255}, 
pages = {106635}, 
year = {2021}, issn = {0045-7949}, 
doi = {https://doi.org/10.1016/j.compstruc.2021.106635}, 
url = {https://www.sciencedirect.com/science/article/pii/S0045794921001577}, 
}

@Article{Rocha2018,
  author    = {Rocha, Felipe Figueredo and Blanco, Pablo Javier and S{\'a}nchez, Pablo Javier and Feij{\'o}o, Ra{\'u}l Antonino},
  title     = {{Multi-scale modelling of arterial tissue: Linking networks of fibres to continua}},
  journal   = {Computer Methods in Applied Mechanics and Engineering},
  year      = {2018},
  volume    = {341},
  pages     = {740--787},
  issn      = {00457825},
  doi       = {10.1016/j.cma.2018.06.031},
  publisher = {Elsevier B.V.},
}


@PhdThesis{FelipeThesis,
  author    = {Rocha, F.F},
  title     = {Multiscale modelling of fibrous materials: from the elastic regime to failure detection in soft tissues},
  school    = {Laborat{\'o}rio Nacional de Computa{\c c}{\~a}o Cient\'{\i}fical},
  year      = {2019},
  publisher = {Laborat{\'o}rio Nacional de Computa{\c c}{\~a}o Cient\'{\i}fica},
}
