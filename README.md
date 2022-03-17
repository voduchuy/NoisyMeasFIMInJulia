# NoisyMeasFIMInJulia

Julia translation of the numerical analyses in the [preprint](https://doi.org/10.1101/2021.05.11.443611):

_Designing single-cell experiments to harvest fluctuation information while rejecting measurement noise_
Huy D. Vo, Brian Munsky
bioRxiv 2021.05.11.443611;

## Setting up the Julia environment
A `Project.toml` file is provided in the root folder to facilitate setting up a Julia environment. 
Simply change working directory to the local folder of the repository. Open a Julia REPL session. Then type `]` to 
switch to package mode. Then enter the following commands
```
pkg> activate
pdg> instantiate
```
The appropriate dependencies (including the pure Julia implementation of the Finite State Projection, [NumCME](https://github.com/voduchuy/NumCME.jl)) 
should be automatically downloaded and installed to the environment.

The user can then execute the scripts provided in each subfolder to produce the numerical outputs, then run the IJulia
notebook to visualize those results.

## Todo
The translation is in progress. The following examples in the preprint are still missing:
- Missing examples related to the 2-state telegraph model:
  - [ ] Noisy cell segmentation FIMs.
  - [ ] Plotting for Binning distortion FIMs.
  - [ ] MLE validation results depicted in Fig.2 in main text.
  - [ ] MLE validation results depicted in Supplemental Fig. S7.
- [ ] All outputs related to toggle-switch model.
- [ ] All outputs related to the MAPK-activated gene expression mode.
