EponaGC
=======

Welcome to the invasive spread modeling program, EponaG(graph)C(control). 
Author: Christopher Strickland (2013)

Epona is a Python 2.7/Cython implementation of the model described in Strickland, Dangelmayr, Shipman (2014) "Modeling the presence probability of invasive plant species with nonlocal dispersal", *Journal of Mathematical Biology*, 69(2), 267--294. The network model, coupling and control in this software is an implementation of the model and optimal control described in Strickland, C. (2013), "The Mathematical Modeling and Analysis of Nonlocal Ecological Invasions and Savanna Population Dynamics", *Colorado State University*, PhD dissertation.

-----------------------------------------------------------------------------------------------------
-----------------------------------------Important!! ------------------------------------------------

This code was created while I was a graduate student, learning Python. I make no claim as to it's
correctness or utility for any problem outside the project it was specifically created for. 
**It is not currently being developed, and you are on your own if you choose to use it for your own project.**

Dependencies
------------

* NumPy
* SciPy
* Matplotlib
* xlrd (for reading excel files)

Installation: compiling the solver (necessary unless running Windows x64)
--------------------

IF YOU DO NOT PLAN ON EDITING THE SOURCE CODE OF THE SOLVER AND YOU ARE RUNNING WINDOWS x64, IGNORE ALL OF THE FOLLOWING. I'VE INCLUDED THE BINARY IN THE win_x64 FOLDER (EponaGSolverC.pyd). JUST COPY IT OVER INTO THE SAME FOLDER AS EponaGC.py, AND LIVE HAPPILY EVER AFTER.

In order to run this program, you will need a compiled version of EponaGSolverC.pyx that matches your operating system, placed in the same folder as EponaGC.py. I've included one for Windows x64 in win_x64, as that is the operating system I'm currently developing on (I know, I know...). Otherwise, assuming Cython is installed on your computer, along with a C compiler that is compatible with Python 2.7, you can do compile using the included setup.py. The command is:  
python setup.py build_ext --inplace

If you are running Windows, and depending on your Python installation, this can be a larger problem than just installing Cython and finding some version of gcc - you must use the C compiler that was originally used to build Python 2.7, and it is no longer supported by Microsoft. Also, old Microsoft installations of this environment do not tend to play well with modern systems. This is where Anaconda comes to your rescue. This Python package not only comes with Cython right out of the box, but also sets up an environment that includes the compiler you need, and makes it just work. Be sure you uninstall/unregister any current Python installations before downloading and installing Anaconda.

How to run simulations with EponaGC
-----------------------------------

Most of the setup for a simulation occurs in the Param.ini file and, in the case of network interactions, the Graph.ini file. You can also apply control on the continuous spatial domain using Cntrl.py. Once these are to your liking, simply run the EponaGC.py script. Press return at the menu, and you will be prompted as to whether or not a transportation network is coupled in the simulation. If yes, information about the network will be pulled from Graph.ini. You will also be given the option to use control as specified in Cntrl.py - Cntrl.py is ignored if n is selected. The simulation will then proceed to run, reporting the current time solved and the current time step, which will be adjusted automatically as the stiffness of the problem allows.

If you would like to begin by reproducing the results of Strickland et al., all the necessary data has been provided in data/data.zip. Just use the default parameters.

### Param.ini and importing data

All non-network information about the current simulation is specified here, including the location of necessary data. An explanation of each of the entries follows.

* `runtime`  Number of years the simulation should be run. Must be an integer.
* `rate`  The underlying population dynamics follow a logistic growth model. This is the growth rate for that model.
* `K`  The carrying capacity, per spatial cell, of the logistic growth model.

* `out_prefix`  What you want the output to be called. EponaGC saves results in .dat cPickle files, the exact format of which will be described later in this document.

* `suitability_data`  The location of an ascii file containing suitability data, space delimited. E.g., maxentmodel08.asc, or data\maxentmodel08.asc. EponaGC was built to run on standard asc output from the MaxEnt suitability model, including a six line header (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value). If the header format is recognized, this information will also be read in by EponaGC.py.

* `initial_data`  The location of a comma delimited (.csv) file specifying initial presence location points for the invader, given in UTM. There should be one column for UTMEasting and one column for UTMNorthing, with each row specifying a different presence point. This file does not have to contain only this data, and may also contain additional presence locations for later years, used for model validation.
* `initdata_range`  Four comma delimited integers that tell EponaGC which rows/columns to read in the initial_data .csv file. In order, they specify: first row to be read, last row to be read, the column with UTMEasting data, the column with UTMNorthing data. It is assumed that row/column numbering begins with 1.

* `pres_data`  Additional presence data can be plotted along with the model for validation. This variable specifies the location of that data, similarly to `initial_data` above. It may be the same file.
* `presdata_col`  Just like the second two values in `initdata_range`, this should be two comma delimited integers specifying the UTMEasting and UTMNorthing columns, in that order.
* `presdata_range`  Comma delimited integers, of which there must be an even number. This list should be thought of as consisting of pairs, each specifying the row range (first row, last row,...) for different sets of presence data to be added to the solution plot at (possibly) different times. This is just like the first two values in `initdata_range`, only there can be multiple "first row, last row" blocks.
* `pres_years`  When to add each block of presence points to the plot. This is a list of integers, and there should be N/2 of them, where N is the number of integers in `presdata_range`. So, 3,5,7 corresponds to "plot the first block of presence locations starting at t=3, the second block at t=5, and the third block at t=7". Only the most recent set is plotted, but note that the initial presence points are always plotted.

* `road_data`  Location of optional ascii file containing raster data for road locations. This is not used in the model - it is for plotting purposes only.

* `negval`  Any presence probability below this number will be assumed to be zero. This prevents your output from filling up with negligible presence probability. It should be a small number, like 0.01.

* `weight`  The type of probability density function to use for the spread kernel. You currently have the option of Normal or Laplace.
* `mu`  Two numbers, comma delimited, giving the mean of the spread kernel in the x direction and y direction, respectively.
* `sig2`  Two numbers, comma delimited, giving the variance of the spread kernel in the x direction and y direction, respectively.
* `corr`  One number, giving the correlation of the spread kernel.

* `correction`  This is by far the most difficult of the parameters in this model - it is the correction value "gamma" given in equation (21) of Strickland et al. "Modeling the presence probability of invasive plant species with nonlocal dispersal." The model can be quite sensitive to this number, and its value is highly dependent on the specifications of the spread kernel and the carrying capacity, K, that was chosen. The best way to determine it is to run basic 1-D model simulations and compare with psuedo-random results from the stochastic process, or check the paper to see if has already been calculated. At this time, EponaGC does not automate the process of determining this value, though it may be added in future development.

You can specify some presence uncertainty around a presence location point in the initial conditions. This is done by a bivariate normal distribution with the following parameters:

* `mu_0` mean in x and y direction, comma delimited.
* `sig2_0` variance in the x and y direction, comma delimited.
* `corr_0` correlation value

If spatial control is to be applied in the model, it is done so via a function specified in Cntrl.py. The following parameters control the effectiveness and capacity of that control.

* `cntrl_K`  maximum total application of control per year
* `cntrl_r`  effectiveness of control against invader
* `cntrl_Lr`  effectiveness of control against seeds, particularly when coupled with network model.

### Graph.ini and specifying a network

All the information about the coupled network is specified here, including parameters for any control on that network (see Strickland, C. dissertation).

* `graphfile`  Location of comma delimited csv file containing:
1. Geographical node locations. One column for UTMEasting, one column for UTMNorthing.
2. Matrix of directed edge weights. Diagonal entries will be ignored and autogenerated, but must not be empty.
3. A column specifying the rate at which individuals leave the network from each node, to come back as susceptible vectors. (mu)
4. A column of probabilities, must sum to 1. When leaving and returning to the network as susceptible vectors, this column represents the probability that they re-enter at each of the nodes. (nu)

The following variables specify how `graphfile` should be read.
* `node_range`  first row, last row, UTMEasting column, UTMNorthing column for geographical node locations.
* `out_rate_range`  first row, column for the rate at which individuals leave (mu)
* `in_prob_range`  first row, column for probabilities governing where they re-enter (nu)
* `graph_range`  first row, first column for the matrix of directed edge weights

The following parameters govern how the network is attached to the spatially continuous model

* `weight`  Probability distribution to be used. Normal and Laplace are supported.
* `mu`  Mean in the x and y direction of this distribution.
* `sig2`  Variance in the x and y direction of this distribution.
* `corr` Correlation value of the distribution

The following parameters turn control on the network on/off, and specify some parameters on it. Control on the network is always the optimal control.

* `cntrl`  1 or 0 (on or off).
* `norm`  integer, l-norm with which to bound amount of control.
* `K`  Upper bound on the norm of the control.

More parameters

* `days`  This deals with the two time-scales between the network and the spatially continuous model. Basically, it should be the length of the growing season in days.
* `decay_rate`  Fraction of seeds that should be lost from the system each year
* `gam_rate`  Rate governing how infectious the spatially continuous model is to network individuals
* `r_rate`  Rate governing how many seeds an infected network individual deposits to the spatially continuous model
* `sprout_rate`  Rate at which seeds become mature invaders.
* `N`  Number of individuals are on the network (S+C). Constant.
* `Lnum`  Allows the user to initialize with some seeds already dropped from the network at node locations.

Output and plotting
-------------------

Output is in cPickle .dat files, one for the spatially continuous result, the other for the dynamics on the network. Because spatial domains are often large, results for the spatially continuous model are only recorded after every integer year. The processed suitability data is also output for the purposes of faster recall, also in a cPickle format, as suitdata.dat. This name can safely be changed after saving - EponaGC gives the option of loading a suitdata file for plotting with a different name.

* The contents of the main .dat file are (in order):
1. model output
2. initial presence data
3. run time
4. x domain info
5. y domain info
6. boundary location (inferred from suitability data)
7. road data
8. non-initial presence data
9. corresponding years for the non-initial presence data
10. locations of the network nodes.

* The contents of the graph .dat file are (in order):
1. susceptibles
2. carriers
3. latent seeds
4. control

The model solution can be plotted right after solving, or from the .dat files. In the latter case, choose either "Plot previous solution" or "Plot previous graph states" as appropriate from the main menu (first thing you see when running EponaGC.py). The first option plots the model in the spatially continuous domain, while the second plots the network dynamics. You will be prompted for a .dat file to load.

In addition to plotting the spatially continuous model at each year, you are given the option to scale the presence probability by suitability and/or plot the suitability separately.
