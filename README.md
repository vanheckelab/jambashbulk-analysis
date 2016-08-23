jambashbulk-analysis
====================
(formerly known as `phd-library`).

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60687.svg)](http://dx.doi.org/10.5281/zenodo.60687)


This is a set of data analysis scripts to work with the output of
the [jambashbulk](https://github.com/vanheckelab/jambashbulk) simulation,
plus some general tools to aid in plotting.

The general overview of this repository is as follows:
 * `packing_tools` contains code to load the output of `jambashbulk`
   and to calculate static properties (contact network, etc) and elastic
   response,
 * `hdf_tools` contains code to write and read simulations in the hdf5 (`.h5`)
   storage format, which is much more efficient in terms of i/o.
 * `upar_uperp` contains code to calculate the full `u_\perp` and `u_\parallel`
   elastic response to a deformation
 * the root directory contains some general tools.

## Context
For background and the physics behind the simulation, please see:
 * Simon Dagois-Bohy, Brian P. Tighe, Johannes Simon, Silke Henkes, and Martin van Hecke.
   _Soft-Sphere Packings at Finite Pressure but Unstable to Shear_.
   [Phys. Rev. Lett. **109**, 095703](http://dx.doi.org/10.1103/PhysRevLett.109.095703),
   [arXiv:1203.3364](http://arxiv.org/abs/1203.3364).

* Merlijn S. van Deen, Johannes Simon, Zorana Zeravcic, Simon Dagois-Bohy, Brian P. Tighe, and Martin van Hecke.
  _Contact changes near jamming_.
  [Phys. Rev. E **90** 020202(R)](http://dx.doi.org/10.1103/PhysRevE.90.020202), [arXiv:1404.3156](http://arxiv.org/abs/1404.3156).

* Merlijn S. van Deen, Brian P. Tighe, and Martin van Hecke.
  _Contact Changes of Sheared Systems: Scaling, Correlations, and Mechanisms_.
  [arXiv:1606.04799](https://arxiv.org/abs/1606.04799)

* Merlijn S. van Deen.
  _Mechanical Response of Foams: Elasticity, Plasticity, and Rearrangements_.
  PhD Thesis, Leiden University, 2016. [hdl:1887/40902](https://openaccess.leidenuniv.nl/handle/1887/40902)

Several data sets created using the simulation code are available via Zenodo:
 * http://dx.doi.org/10.5281/zenodo.59216 (stable packings)
 * http://dx.doi.org/10.5281/zenodo.59217 (sheared packings)

## Contents

### General tools
 * [cdf.py](cdf.py) contains tooling to plot CDFs,
 * [smart_annotate.py](smart_annotate.py) contains code to click-and-mark
   data points in matplotlib graphs,
 * [util.py](util.py) contains helper functions to find external hard drives
   and to consistently use the same markers/colors for plotting.

### Packing tools
 * [load_packing.py](packing_tools/load_packing.py) contains logic to parse the
   `jambashbulk` packing data format (`N16~P1e-3~0001.txt`, `particlesN16~P1e-3~SR003~step007~0001.txt`)
   * [parser](packing_tools/parser) is a C-based parser to speed up this process.
 * [load_log.py](packing_tools/load_log.py) parses `logN16~P1e-3~0001.txt`.
 * [make_shear_graphs.py](packing_tools/make_shear_graphs.py) plots packings,
 * [V_harm.py](packing_tools/V_harm.py) calculates the contact network, Hessian, elastic moduli, etc.

### HDF tools
 * [pytables_test.py](hdf_tools/pytables_test.py) loads static packings into `N16~P1e-3_tables.h5` files
 * [pytables_import_shear.py](hdf_tools/pytables_import_shear.py) loads contact change data into `N16~P1e-3_shear.h5` files
 * [generate_cache.py](hdf_tools/generate_cache.py) generates an hdf5 file with information on the first contact change
   in each simulation
 * [dodo.py](dodo.py) (in the root) automates these three steps for a large number of packings.
 * [pytables_tools.py](hdf_tools/pytables_tools.py) has the `read_packing` function to read packings from the HDF format
 * [getcc.py](hdf_tools/getcc.py) contains code to parse individual contact changes from the shear simulation data


### upar/uperp
 * [hpc.py](upar_uperp/hpc.py) contains the `HessianPackingCalculator` calculates upar and uperp from
   linear response, and calculates the strain at which the first contact change happens
   from this information.
