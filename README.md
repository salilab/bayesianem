# bayesianem, an IMP module for Bayesian multi-scale modeling of macromolecular structures based on cryo-electron microscopy density maps

[![Build Status](https://github.com/salilab/bayesianem/actions/workflows/build.yml/badge.svg)](https://github.com/salilab/bayesianem/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/salilab/bayesianem/branch/main/graph/badge.svg)](https://codecov.io/gh/salilab/bayesianem)

## Usage:

First, import the module:

```
import IMP.bayesianem
import IMP.bayesianem.restraint
```

Then, select the model densities that you want to assign to the EM:

```
sel = IMP.atom.Selection(hierarchy=hier, representation_type=IMP.atom.DENSITIES)
densities=sel.get_selected_particles()
```

Finally, create the restraint (replace gmm.txt by your target gmm), and add it to the model.

```
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                 target_fn=‘gmm.txt',
                                                 scale_target_to_mass=True,
                                                 slope=0.01,
                                                 target_radii_scale=3.0,
                                                 target_is_rigid_body=False)
gem.add_to_model()
gem.set_label("Total")
output_objects.append(gem)
```

## Info

_Author(s)_: Samuel Hanot, Massimiliano Bonomi, Riccardo Pellarin

_Maintainer_: Riccardo Pellarin

_License_: LGPL

_Publications_:
 - M. Bonomi, S. Hanot, C. H. Greenberg, A. Sali, M. Nilges, M. Vendruscolo, R. Pellarin. ["Bayesian Weighing of Electron Cryo-Microscopy Data for Integrative Structural Modeling", Structure 27, 2019](https://pubmed.ncbi.nlm.nih.gov/30393052/).
