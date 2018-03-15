# SBD-iPALM: Sparse blind deconvolution using iPALM

**SBD-iPALM** is a MATLAB package for *sparse blind deconvolution* (SBD) using the iPALM method, motivated by studies in blind deconvolution as a nonconvex optimization problem, and by applications in Scanning Tunneling Microscopy (STM).

The [iPALM method](https://arxiv.org/abs/1702.02505) is an accelerated first-order method with optimal convergence guarantees to a stationary point when solving nonconvex optimization problems. We combine this with the [l1-reweighting method](https://arxiv.org/abs/0711.1612) to produce sparse, robust and interpretable activation maps.

Please contact [yl3027@columbia.edu](yl3027@columbia.edu) for any requests / feedback. We thank Sky Chueng and Abhay Pasupathy for the simulated kernel data `Data_N_50_nDef_1_thop_-0.200.mat` used in `genconvdata.m` and `example.m`.


## Last update for v2
See [UPDATES.md](UPDATES.md) for complete list.

**2018-03-16**
- Some modifications to make SBD-iPALM more useful as a standalone package

    * `initpkg.m`: adds all utilities / scripts to the MATLAB path

    * Some default iteration loops and updates are provided as seperate scripts

        * `reweight_loop.m`:  default looping of iterations with reweighting

        * `text_update.m`:  default update text print

        * `cdl_update.m`:  default SBD / CDL centering, plots, update text

    * Fixed centering for general CDL problems.



# Expected additions

- An ADMM iterator may be added in the future.

- Reimplement wrappers of SBD, CDL, etc. as classes? To add various associated methods such as synthesis.