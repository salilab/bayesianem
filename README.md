# bayesianem, an IMP module for Bayesian multi-scale modeling of macromolecular structures based on cryo-electron microscopy density maps

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
Bayesian multi-scale modeling of macromolecular structures based on cryo-electron microscopy density maps
Samuel Hanot, Massimiliano Bonomi, Charles H Greenberg, Andrej Sali, Michael Nilges, Michele Vendruscolo, Riccardo Pellarin
bioRxiv 113951; doi: https://doi.org/10.1101/113951 (submitted)
