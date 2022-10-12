# Nexis
Repository Owner: [Justin Torok](http://github.com/justin-torok)
Email: justin.torok@ucsf.edu

Project Lead: Chaitali Anand
Email: chaitali.anand@ucsf.edu

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
- `stdNDM_mouse.m`: The core function used to generate Nexis:global results (refer to [Anand, *et al.*, 2022](https://www.biorxiv.org/content/10.1101/2021.03.22.436470v1) for full model details). This function solves the parameter inference problem and saves final model outputs with summary statistics in a MATLAB struct object. All inputs are optional and specified as keyword arguments using `inputParser`. 
    - ***Inputs***:
        - `study`: String indicating which dataset to use (default = 'IbaHippInj')
        - `costfun`: String indicating which cost function to use, which is passed to `objfun_eNDM_general_dir_costopts.m` (default = 'LinR')
        - `solvetype`: String indicating whether to use an analytical solution to the Nexis:global model or a numerical implentation (default = `analytic`)
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
