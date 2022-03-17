If this is the first time the user tries to re-run the analyses, please make a folder called `results`
in this subfolder
```zsh
mkdir results
```
We will store all numerical outputs (as `.jld2` files) in this subfolder.

## Script execution
We assume that a local Julia environment has been appropriately instantiated with all dependencies installed (see Readme file in the root folder). The very first script to execute
```
julia --project=<path_to_root_folder> 0_solve_sensitivity_cme.jl
```
These generate the model-predicted mRNA distributions and their partial derivatives w.r.t model parmeters. Subsequent scripts simulate how parameter uncertainties are distorted under different measurement errors.

## Note on the order of script execution 
The scripts are named in a hiearachical way. 
The script `2_***.jl` should not be executed until `i_***.jl` for 
i in [0,1] have all been executed since it may depend on the outputs at those preceding scripts. 
Scripts with the same numerical prefix can be executed in 
any order relative to each other.

The user is advised that the script `1_flowcyt_FIMs.jl` will take about 30min to execute.

## Result visualization
Once all scripts have been executed to populate the `results` subfolder, the user can run `bursting_gene_plots.ipynb` to 
plot them. These plots have similar layouts to the original Matplotlib-generated plots in the preprint. However, due to 
differences between Matplotlib and Julia Plots, we cannot reproduce the figures 100% pixel by pixel.
