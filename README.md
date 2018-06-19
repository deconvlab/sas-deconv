# SBD-iPALM: Sparse blind deconvolution using iPALM

**SBD-iPALM** is a MATLAB package for *sparse blind deconvolution* (SBD) using the iPALM method, motivated by studies in blind deconvolution as a nonconvex optimization problem, and by applications in Scanning Tunneling Microscopy (STM).

The [iPALM method](https://arxiv.org/abs/1702.02505) is an accelerated first-order method with optimal convergence guarantees to a stationary point when solving nonconvex optimization problems. We combine this with the [l1-reweighting method](https://arxiv.org/abs/0711.1612) to produce sparse, robust and interpretable activation maps.

Please contact [Yenson Lau](yl3027@columbia.edu) for any requests / feedback. We thank [Sky Chueng](physics.columbia.edu/people/profile/469), [John Shin](physics.columbia.edu/people/profile/938) and [Abhay Pasupathy](physics.columbia.edu/people/profile/428) for the simulated kernel data `Data_N_50_nDef_1_thop_-0.200.mat` used in `genconvdata.m` and `example.m`.


## Most recent update for v2
**2018-06-20**
- Added `PDRegularizer` for regularizers of the form g(Dx)

    * comes with a Chambolle-Pock primal-dual method to solve for prox updates 

- A difference operator `imgdiff` for image deblurring: plug into D above. Improvements in progress.

- Numerous minor adjustments to clean up the code, which will continue to roll in the next few updates.

See [UPDATES.md](UPDATES.md) for complete list.


## Expected additions

- An ADMM iterator may be added in the future.

- Reimplement wrappers of SBD, CDL, etc. as classes? To add various associated methods such as synthesis.