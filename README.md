# SBD-iPALM: Sparse blind deconvolution using iPALM

**SBD-iPALM** is a MATLAB package for *sparse blind deconvolution* (SBD) using the iPALM method, motivated by studies in blind deconvolution as a nonconvex optimization problem, and by applications in Scanning Tunneling Microscopy (STM).

The [iPALM method](https://arxiv.org/abs/1702.02505) is an accelerated first-order method with optimal convergence guarantees to a stationary point when solving nonconvex optimization problems. We combine this with the [l1-reweighting method](https://arxiv.org/abs/0711.1612) to produce sparse, robust and interpretable activation maps.

Please contact [yl3027@columbia.edu](yl3027@columbia.edu) for any requests / feedback. We thank Sky Chueng and Abhay Pasupathy for the simulated kernel data `Data_N_50_nDef_1_thop_-0.200.mat` used in `genconvdata.m` and `example.m`.


## Updates for v2

**2018-03-13 (First batch of GitHub commits `GitHub1`)**:

- Improved usability + commenting for core scripts:
    * `ipalm.m`, `mksbd.m`,  `mkcdl.m`, `genconvdata.m`

- `example.m` updated. Comments + centering fixed.


**2018-03-08**:

- `mkcdl.m` updated to work with multiple samples. Centering needs to be updated and commenting still needed.


**2018-03-05**:

- Generic iterators for solving various bilinear problems on the sphere -- e.g. `ipalm` -- now power problem instances of interest -- e.g. SBD, CDL.

- SBD and CDL now implemented as wrapper functions for constructing specific iterator instances.

- The `huber` loss is now rewritten to acommodate for the iterator structure. Reweighting is included.

- All iterating, reweighting and display functions that go on top of the inner loop can now be written seperately, e.g. this means one can run iterations until desired and not lose information. Changes to each problem setting can be made easily.


# Expected updates

- An ADMM iterator may be added in the future.

- Reimplement wrappers of SBD, CDL, etc. as classes? To add various associated methods such as synthesis.