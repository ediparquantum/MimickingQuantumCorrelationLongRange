# Mimicking quantum correlation of a long-range Hamiltonian by finite-range interactions

The associated code to the paper [Mimicking quantum correlation of a long-range Hamiltonian by finite-range interactions](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.052425) written by Leela Ganesh Chandra Lakkaraju, Srijon Ghosh, Debasis Sadhukhan, and Aditi Sen(De) and was published 22 November 2022 by PHYSICAL REVIEW A. 

## Abstract

The quantum long-range extended Ising model possesses several striking features which cannot be observed in the corresponding short-range model. We report that the pattern obtained from the entanglement between any two arbitrary sites of the long-range model can be mimicked by the model having a finite range of interactions provided the interaction strength is moderate. On the other hand, we illustrate that when the interactions are strong, the entanglement distribution in the long-range model does not match the class of a model with a few interactions. We also show that the monogamy score of entanglement is in good agreement with the behavior of pairwise entanglement. Specifically, it saturates when the entanglement in the finite-range Hamiltonian behaves similarly to the long-range model, while it decays algebraically otherwise.

## Code
Code relies on the [Quantum Information and Computation library (QIClib)](https://titaschanda.github.io/QIClib/). 

### Running the code

1.  Run `src/execute_code_submit0.sh` in a terminal (Linux) so that Armadillo is findable (Note it is possible this script is relevant to a local machine and will need to be upated based on user defined parameters).
2. Run the file `src/deba_Z_h.cpp`

### Data
Data output is found at `data/xZ_g_norm_1_n512.dat`


### QIClib

The QIClib library has the following benefits.

+ Mordern C++11 library suited for general purpose quantum computing. 
+ It supports cross platform usage (Linux, Windows and Mac OS X). 
+ It is a header only template library.
+ It uses Armadillo (developed by Conrad Sanderson et al., Data61, Australia) for highly efficient linear algebra calculations, and if available, the NLopt nonlinear optimization library for certain features.
    + [Armadillo home page](https://arma.sourceforge.net/)

To see QIClib's GitHub goto: [https://github.com/titaschanda/QIClib](https://github.com/titaschanda/QIClib)


## Citation to paper

```
@article{PhysRevA.106.052425,
  title = {Mimicking quantum correlation of a long-range Hamiltonian by finite-range interactions},
  author = {Lakkaraju, Leela Ganesh Chandra and Ghosh, Srijon and Sadhukhan, Debasis and Sen(De), Aditi},
  journal = {Phys. Rev. A},
  volume = {106},
  issue = {5},
  pages = {052425},
  numpages = {13},
  year = {2022},
  month = {Nov},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevA.106.052425},
  url = {https://link.aps.org/doi/10.1103/PhysRevA.106.052425}
}
```

## Accessing the paper

A [pdf version as found on arXiv](https://arxiv.org/abs/2206.09199) is uploaded to the repository.

+ `docs/arXiv-Can a finite range Hamiltonian mimic quantum correlation of a long-range Hamiltonian.pdf`



