# Nexis
Repository Owner: [Justin Torok](http://github.com/justin-torok) (Email: justin.torok@ucsf.edu)

Project Lead: Chaitali Anand (Email: chaitali.anand@ucsf.edu)

The following code was developed for running the Nexopathy in silico (NexIS) model described in Anand, *et al.*, 2022 (preprint [here](https://www.biorxiv.org/content/10.1101/2021.03.22.436470v1)), along with all of the auxilliary functions required for plotting the outputs as shown in the manuscript. 

We note that the naming conventions of the models have changed, but these changes are not reflected in the code at this time. The NexIS:global and NexIS:microglia results described in the manuscript were run using the `stdNDM_mouse.m` and `eNDM_mouse.m` functions, respectively (these are described in more detail below). 

## 1. Setup
All code is written in MATLAB and has been tested in versions 2020a and 2022a. However, the developers do not anticipate difficulties with using other MATLAB versions for any of the functions contained within this package. There are no toolboxes required to run the model, but the Statistics & Machine Learning Toolbox is required for the `corr.m` function. All necessary data dependencies are within the raw_data_mouse folder, which contains the following data germane to the present study:
- Mouse tauopathy data, which were obtained from [Kaufman, *et al.*, 2017](https://pubmed.ncbi.nlm.nih.gov/28587664/), 
- Regional gene expression data, which were obtained from the [Allen Gene Expression Atlas](https://mouse.brain-map.org/), 
- Mouse mesoscale connectome, which was obtained from [Oh, *et al.*, 2014](https://www.nature.com/articles/nature13186). 

## 2. Files
Below is a short description of each of the code files contained in the **Nexis** repository, grouped by general functionality. Scripts that are also functions have their inputs and outputs described, with required inputs in boldface text and optional inputs with their default setting in parentheses.

### Running the Model
- `stdNDM_mouse.m`: The core function used to generate Nexis:global results (refer to [Anand, *et al.*, 2022](https://www.biorxiv.org/content/10.1101/2021.03.22.436470v1) for full model details). This function solves a parameter inference problem on [gamma, alpha, beta, s] and saves final model outputs for the optimized values of these parameters, with summary statistics, in a MATLAB struct object. All inputs are optional and specified as keyword arguments using `inputParser`. 
    - ***Inputs***:
        - `study`: String indicating which dataset to use (default = 'IbaHippInj')
        - `costfun`: String indicating which cost function to use, which is passed to `objfun_eNDM_general_dir_costopts.m` (default = 'LinR')
        - `solvetype`: String indicating whether to use an analytical solution to the NexIS:global model or a numerical implentation (default = `analytic`)
        - `volcorrect`: Binary flag indicating whether to add a volume correction to the model. In practice, this does not affect model results very much (default = 0)
        - `exclseed_costfun`: Binary flag indicating whether or not to exclude the seed regions from the cost function calculation (default = 0)
        - `excltpts_costfun`: List of time point indices to exclude from the cost function calculation (default = [])
        - `normtype`: String indicating which normalization, if any, to use on the data (default = 'sum' - **we recommend changing this to 'none'**)
        - `w_dir`: Binary flag indicating whether or not to use the directional connectome. If 0, the connectome will be symmetrized prior to running the model and the parameter 's' will be fixed at 0.5 (default = 0)
        - `param_init`: Initial parameter values to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model. Unless otherwise specified, the initial value of gamma is determined heuristically in the function itself (default = [NaN,0.5,1,0.5])
        - `ub`: Upper bounds to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model (default = [Inf,Inf,Inf,1])
        - `lb`: Lower bounds to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model (default = zeros(1,4))
        - `algo`: String specifying the `fmincon` algorithm (default = 'sqp')
        - `opttol`: Value of the optimality tolerance in `fmincon` (default = 1e-8)
        - `fxntol`: Value of the function tolerance in `fmincon` (default = 1e-8)
        - `steptol`: Value of the optimality tolerance in `fmincon` (default = 1e-12)
        - `maxeval`: Maximum number of function evaluations in `fmincon` (default = 10000)
        - `bootstrapping`: Binary flag indicating whether the parameters should be fit using bootstrapping on regions (default = 0)
        - `resample_rate`: If `bootstrapping` = 1, the fraction of regions to use to fit the model per iteration (default = 0.8)
        - `niters`: If `bootstrapping` = 1, the number of bootstrapping iterations (default = 100)
        - `verbose`: Binary flag indicating whether to use verbose in-line outputs (default = 0)
        - `fmindisplay`: Binary flag indicating whether to use verbose `fmincon` displays (default = 0)
         - `flowthresh`: Percentile at which to threshold the calculated flows using `FlowCalculator.m` (default = 99.93)    
    - ***Outputs***:
        - `outputs`: A MATLAB struct containing inputs, outputs, and summary statistics of the final model. 
- `eNDM_mouse.m`: The core function used to generate NexIS:Trem2 results, among others (refer to [Anand, *et al.*, 2022](https://www.biorxiv.org/content/10.1101/2021.03.22.436470v1) for full model details). This function solves a parameter inference problem on [gamma, alpha, beta, s] and saves final model outputs for the optimized values of these parameters, with summary statistics, in a MATLAB struct object. All inputs are optional and specified as keyword arguments using `inputParser`. This function uses a hierarchical parameter fitting approach, where requires a run of `stdNDM_mouse.m` is used to first fit the global parameters [gamma, alpha, beta, s] and then uses those outputs to define the initial values and bounds of these parameters to pass to `fmincon` when fitting for the new parameters [b, p]. For iterating over genes or cell types, we highly recommend running `stdNDM_mouse.m` first and then passing its output struct to `eNDM_mouse.m` to save computational time, but the function will call `stdNDM_mouse.m` internally if this struct is not provided.  

    - ***Inputs***:
        - `study`: String indicating which dataset to use (default = 'IbaHippInj')
        - `costfun`: String indicating which cost function to use, which is passed to `objfun_eNDM_general_dir_costopts.m` (default = 'LinR')
        - `solvetype`: String indicating whether to use an analytical solution to the NexIS model or a numerical implentation (default = `analytic`)
        - `volcorrect`: Binary flag indicating whether to add a volume correction to the model. In practice, this does not affect model results very much (default = 0)
        - `exclseed_costfun`: Binary flag indicating whether or not to exclude the seed regions from the cost function calculation (default = 0)
        - `excltpts_costfun`: List of time point indices to exclude from the cost function calculation (default = [])
        - `normtype`: String indicating which normalization, if any, to use on the data (default = 'sum' - **we recommend changing this to 'none'**)
        - `w_dir`: Binary flag indicating whether or not to use the directional connectome. If 0, the connectome will be symmetrized prior to running the model and the parameter 's' will be fixed at 0.5 (default = 0)
        - `param_init`: Initial parameter values to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model. Unless otherwise specified, the initial value of gamma is determined heuristically in the function itself (default = [NaN,0.5,1,0.5])
        - `ub`: Upper bounds to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model (default = [Inf,Inf,Inf,1])
        - `lb`: Lower bounds to pass to `fmincon`. In order, these correspond to gamma, alpha, beta, and s in the model (default = zeros(1,4))
        - `algo`: String specifying the `fmincon` algorithm (default = 'sqp')
        - `opttol`: Value of the optimality tolerance in `fmincon` (default = 1e-8)
        - `fxntol`: Value of the function tolerance in `fmincon` (default = 1e-8)
        - `steptol`: Value of the optimality tolerance in `fmincon` (default = 1e-12)
        - `maxeval`: Maximum number of function evaluations in `fmincon` (default = 10000)
        - `bootstrapping`: Binary flag indicating whether the parameters should be fit using bootstrapping on regions (default = 0)
        - `resample_rate`: If `bootstrapping` = 1, the fraction of regions to use to fit the model per iteration (default = 0.8)
        - `niters`: If `bootstrapping` = 1, the number of bootstrapping iterations (default = 100)
        - `verbose`: Binary flag indicating whether to use verbose in-line outputs (default = 0)
        - `fmindisplay`: Binary flag indicating whether to use verbose `fmincon` displays (default = 0)
        - `outputs_ndm`: Struct object from `stdNDM_mouse.m`. If not provided, `stdNDM_mouse.m` will called to get estimates of the initial values of the global parameters (default = [])
        - `bounds_type_endm`: String indicating which type of upper and lower limits to pass to `fmincon` based on the bootstrapped distributions of the global parameters from `stdNDM_mouse.m`, if bootstrapping is run. 'CI_X' uses the Xth confidence interval, while 'old' uses +/- 30% of the mean parameter values (default = 'old') 
        - `bootstrapping_endm`: Binary flag indicating whether the parameters should be fit using bootstrapping on regions (default = 0)
        - `resample_rate_endm`: If `bootstrapping` = 1, the fraction of regions to use to fit the model per iteration (default = 0.8)
        - `niters_endm`: If `bootstrapping` = 1, the number of bootstrapping iterations (default = 100)
        - `verbose_endm`: Binary flag indicating whether to use verbose in-line outputs (default = 0)
        - `fmindisplay_endm`: Binary flag indicating whether to use verbose `fmincon` displays (default = 0)
        - `datatype_endm`: String identifying whether to use cell types or genes to add to the NexIS. The regional gene expression data come from the AGEA and the regional cell-type density data from from [Mezias, *et al.*, 2022](https://www.pnas.org/doi/10.1073/pnas.2111786119) (default = 'gene')
        - `datalist_endm`: Genes or cell types to add to the model. This can be specified as a list of indices, if these are known, or as a cell array of names. **Note: if using names, a cell array must be used even if there is only a single name in the list** (default = 3578 (AGEA index for Trem2))
        - `datapca_endm`: Binary flag indicating whether to use the first PC scores of a list of 2 or more genes or cell types (default = 0) 
        - `flowthresh`: Percentile at which to threshold the calculated flows using `FlowCalculator.m` (default = 99.93)    
    - ***Outputs***:
        - `outputs`: A MATLAB struct containing inputs, outputs, and summary statistics of the final model. 

- `FlowCalculator.m`: Function that calculates the flows on the connectome graph over time, which use a first-order approximation for the time derivatives (refer to [Anand, *et al.*, 2022](https://www.biorxiv.org/content/10.1101/2021.03.22.436470v1) for full details). **This function is called internally within `stdNDM_mouse.m` and `eNDM_mouse.m` and generally does not need to be called outside of these functions.** 

    - ***Inputs***:
        - `X_`: Inferred pathology values over time as a regions-by-time-points matrix.
        - `C_`: Mouse mesoscale connectome (regions-by-regions)
        - `beta_`: Rate parameter for spread on the connectome, which is obtained from NexIS
        - `isndm`: Binary flag indicating whether to use NexIS:global or NexIS:extended 
        - `U_`: Regions-by-types matrix of gene expression or cell-type density, necessary for NexIS:extended
        - `b_`: NexIS:extended rate parameter(s) 
    - ***Outputs***:
        - `flow_mat`: Regions-by-regions-by-time-points array of calculated flow values

### Interpreting the outputs
- `Output2Table.m`: Produces a MATLAB table of key information from the outputs struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`. Optionally writes to a .csv file.

    - ***Inputs***:
        - `outputs`: MATLAB struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`
        - `writeout`: Optional binary flag indicating whether or not to write to a .csv (default = 0)
        - `filename`: Custom file name string for the table written to a .csv (default contains the date, the study name, and whether it is NexIS:global or NexIS:extended)
        - `filepath`: Directory (string) to which the .csv should be written (default = cd)
    - ***Outputs***:
        - `summarytable`: MATLAB table object that unpacks key information from outputs

- `BootstrappingPlotter.m`: Creates boxplots of the parameter distributions from the outputs struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`. 

    - ***Inputs***:
        - `outputs`: MATLAB struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`
        - `savenclose`: Binary flag indicating whether to write the figure to an uncompressed .tiff file and then close the figure handle (default = 0)

- `CorrelationPlotter.m`: Creates correlation plots per timepoint of modeled pathology vs. expected pathology and, if bootstrapping was used, boxplots of the R<sup>2</sup> distributions across iterations. Uses the outputs struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`. 

    - ***Inputs***:
        - `outputs`: MATLAB struct object from `stdNDM_mouse.m` or `eNDM_mouse.m`
        - `savenclose`: Binary flag indicating whether to write the figure to an uncompressed .tiff file and then close the figure handle (default = 0)

- `CorrelationPlotter_single.m`: Creates correlation plots per timepoint of modeled pathology vs. expected pathology. This function has more inherent control over the way the scatterplots are rendered, but requires some unpacking of the outputs struct from `stdNDM_mouse.m` or `eNDM_mouse.m`. (Will likely be subsumed into `CorrelationPlotter.m` at some point)

    - ***Inputs***:
        - `predicted`: Regions-by-time-points matrix of predicted pathology values
        - `data`: Regions-by-time-points matrix of empirical pathology values
        - `tpts`: 1-by-time-points vector of time points
        - `study`: String handle indicating the tauopathy study used
        - `sv_types`: Cell array of genes or cell types used in fitting the NexIS model (default = {})
        - `geneorct`: String indicating whether genes or cell types were used in fitting the model (default = `genes`)
        - `ispca`: Binary flag indicating whether PCA was used on the gene or cell-type inputs to the model (see `eNDM_mouse.m`; default = 0)
        - `wdir`: Binary flag indicating whether the directionality parameter, s, was fit in NexIS (see `stdNDM_mouse.m` or `eNDM_mouse.m`; default = 0)
        - `color`: 1-by-3 vector of RGB values for the data points in the scatterplot(default = [1 0 0])
        - `shape`: Valid string indicating the data point shape for `scatter.m` to use (default = 'o')
        - `savenclose`: Binary flag indicating whether to write the figure to an uncompressed .tiff file and then close the figure handle (default = 0)