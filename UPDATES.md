# Updates for v2

**2018-06-15**
- Final version of SBD-iPALM v2: starting v3, which will simplify `iPALM.m`

**2018-03-16**
- Some modifications to make SBD-iPALM more useful as a standalone package

    * `initpkg.m`: adds all utilities / scripts to the MATLAB path

    * Some default iteration loops and updates are provided as seperate scripts

        * `reweight_loop.m`:  default looping of iterations with reweighting

        * `text_update.m`:  default update text print

        * `cdl_update.m`:  default SBD / CDL centering, plots, update text

    * Fixed centering for general CDL problems.


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