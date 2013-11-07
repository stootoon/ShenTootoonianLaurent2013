# Introduction

This repository contains the code used to produce the figures in the article "*Encoding of Mixtures in a Simple Olfactory System*" by K. Shen, S. Tootoonian and G. Laurent, published in the journal *Neuron* in November 2013. The associated data, consisting of the raw spike times of the cells analyzed in the paper as well as the processed data used to produce the figures is also provided. The code was written almost entirely in MATLAB 2009b for Linux (Ubuntu) and tested on Linux (Ubuntu) and Mac OS X 10.6 (Snow Leopard). It should also run on Windows with minor modification.

# Installation
To install the code and data,

1. Download the code repository from [this](https://dl.dropboxusercontent.com/u/5517909/ShenTootoonianLaurent2013Code_v0.9.zip) link and unpack it. The top-level of the unpacked repository should contain a `code` folder, a few `.m` files and `README.md`. This will be the **installation root**.
2. Download the data from [this](https://dl.dropboxusercontent.com/u/5517909/ShenTootoonianLaurent2013Data.zip) link. It's ~800 MB, so it might take a while.
3. Unpack the data file to yield a `data` folder.
4. Move this folder into the **installation root**. It should now sit at the same level as the `code` folder.
5. Start MATLAB, cd to the **installation root**, and run `install();` This script will set the root directory of the package and compile a mex file. Note that if you start MATLAB in the **installation root** before running `install()`, MATLAB will run `startup` and will complain about INSTALL_ROOT. This warning can be safely ignored, and will disappear once `install()` has been run once.

6. Run `startup();` This will add some code directories to the path.

You can test the installation by checking the panels available for Figure 2, and plotting the raster in Figure 2A-D:

    	>> cd code/figures/figure2;
		>> MakeFigures;
		Panels available for Figure 2: A B C D E F SA SB SC SD SE SF 

		USAGE:
		Plot one panel e.g. 'A':			MakeFigures('A'); 
		Plot several panels e.g. 'A','SB':	MakeFigures({'A','SB'}); 
		Plot all panels:					MakeFigures('all');
		>> MakeFigures('A');

After a few seconds, a figure should popup with the rasters in Figure 2A-D. Success!

# Folder Contents
## Code
The code tree is rooted at `code`. It contains three folders:

1. `util`: Contains general utility functions.
2. `common`: Contains commonly used functions for analysis and plotting of the data.
3. `figures`: Contains directories containing code for figures 2-8. 
4. `figures/figures[2-8]`: Each directory contains 
	1. `MakeFigures.m`: Used to plot the panels for the figure.
	2. Zero or more subdirectories containing code for computing and plotting specific panels. 
	3. `ProcessData.m`: Contained by each subdirectory, and used to recompute all the of the data for the associated panel. 
 	 
## Data
The data tree is rooted at `data`. It also contains three folders:

1. `spt`: Contains the raw spike times data. 
2. `odors`: Contains data about the odors used in the complex mixtures experiments.
3. `proc/figures`: Contains directories for figures 2-8, each containing processed data used to plot the associated figure.
4. `proc/figures/figure[2-8]`: Contains two subdirectories:
	1. `originalData`: Post-processed data as provided in this package and used to plot the panels in the paper.
	2. `recomputedData`: Post-processed data recomputed by the the `ProcessData` functions. 

# Datasets
Spike times from two sets of experiments were analyzed in the paper and have been provided here.
## Spike Times
The spike times provided in `data/spt` are matrices in 'toc' (trials-odors-cells) format: Each column contains the spikes times of the response of one cell in one trial to one odor. The columns are arranged by iterating over trials, then odors, then cells. For example, if there are T trials, R odors, and C cells, 

* Columns 1 to T are the response of cell 1 to odor 1 in consecutive trials, 
* Columns T+1 to 2T are its response to odor 2, 
* Columns (R-1)T + 1 to RT are its response to odor R, 
* Column RT+1 is the response of cell 2 to the first trial of odor 1,  

and so on. Equivalently, if A is a 4-D matrix such that A(n,i,j,k) is the time of the n'th spike in the response in trial i of odor j of cell k, then 

		Atoc = reshape(A, [], T * R *C);
		A	 = reshape(Atoc, [], T, R, C);

The columns of the the toc matrices contain the actual spike times. If a cell produced fewer spikes than the length of the column, the remainder of the column is set to zero. This could in principle be problematic if cells also tend to produce spikes at exactly t = 0, but this is very rare because odor onset is at t = 2.0 seconds and the baseline period is long (> 10 seconds). The function `ConvertSpikeTimesFromSparseToFull` converts the sparse matrices to full by setting the 0 elements to -Inf.

## Binary Mixtures
In this set of experiments, the responses of 168 PNs to 10 trials of 27 mixtures of octanol and citral were recorded. Odor onset was at t = 2.0 seconds and the odor duration was 0.3 seconds. The function `GetBinaryMixturePairedConcentrations` returns the list of (octanol, citral) mixtures used, in the order they are stored in the spike times matrix. To get the PN spike times for this experiment, use:

		>> pnSpt = LoadTocSpikeTimes('rawpn_binary_mixtures');

## Complex Mixtures
In this set of experiments, the responses of 174 PNs and 209 KCs to 7 trials to 8 monomolecular odors and their mixtures was recorded, for a total of 44 odors tested. Odor onset was at t = 2.0 seconds and the odor duration was 0.5 seconds. The function `GetOdorsList` returns the list of mixtures used in the order they are stored in the spike times matrices. To get the spike times for this experiment, use:

		>> pnSpt = LoadTocSpikeTimes('rawpns');
		>> kcSpt = LoadTocSpikeTimes('rawkcs');

Late in the analysis it was discovered that due to a spike-sorting error, one of the PNs in the data set (which actually contains 175 cells) was a duplicate. The affected analyses were rerun with this cell removed, hence the occasional '_noduplicates' suffix of some of the files and folders.

# Basic Usage
The aim of this package is to allow the user to reproduce the figures in the paper. Since most of the figures plot the results of computations performed on the raw spike times, the package also allows the user to reproduce these computations.
## Startup
Always make sure to run the `startup` script in the **installation root** before running the rest of the code, as this adds necessary utility code folders to the path.
## Browsing Spike Times
The `SpikeTimesBrowser` has been provided to allow easy browsing of the spike times provided. Some typical usage examples are shown below; see the function documentation for more details.

        >> SpikeTimesBrowser();

This will plot the raster showing the response of PN 1 in the binary mixtures experiments, in the format of Figure 2A-D. The left and right arrowkeys will show the rasters for previous and subsequent PNs, respectively.

        >> SpikeTimesBrowser('startingCell', 128, 'showRecs', true);

This will reproduce Figure 2A-D, by plotting the raster for PN 128, and showing the reconstructions of its PSTH using linear summation of the component responses.

        >> SpikeTimesBrowser('experiment', 'ComplexMixtures');
        
This will plot the raster for the first PN in the complex mixtures experiments, in the format of Figure 6A-B.

        >> SpikeTimesBrowser('experiment', 'ComplexMixtures', 'cells', 'KCs');

This will plot the raster for the first KC in the complex mixtures experiments.

        >> SpikeTimesBrowser('experiment', 'ComplexMixtures', 'cells', 'KCs', 'startingCell', 155, 'startTime', -0.3, 'endTime', 2.2);

This will reproduce the raster in Figure 6B.

## Reproducing Figures
To reproduce a figure, first `cd` to its code directory. For example, to plot some panels from Figure 3,

		>> cd code/figures/figure3;

You can list all the panels available for this figure by calling `MakeFigures` without any arguments:

		>> MakeFigures;
		Panels available for Figure 3: A B C D E F G H I SA SB SC SD 

		USAGE:
		Plot one panel e.g. 'A':			MakeFigures('A'); 
		Plot several panels e.g. 'A','SB':	MakeFigures({'A','SB'}); 
		Plot all panels:					MakeFigures('all'); 

We can then plot one of the panels, say 3G, the PAF w.r.t citral time course:

		>> MakeFigures('G');

Supplementary panels are prefixed with S. So to plot S3A, the PMF time course:

		>> MakeFigures('SA');

We can plot several panels at once by providing them as a cell array:

		>> MakeFigures({'SA','G'});

Or we can just plot all available panels:

		>> MakeFigures('all');

## Reproducing Computations
To reproduce the computations for a figure panel, first determine the associated subdirectory by using `GetFolderForPanel`. For example, for panel 3G:

		>> GetFolderForPanel('3G');

Then `cd` to that directory and run `ProcessData`. This will perform all the computations required for the panel, writing the results to the associated `recomputedData` folder in the data tree. The results of the computation can then be used to recreate the plot:

		>> cd code/figures/figure3;
		>> MakeFigures('G','dataDir','recomputedData'); 


# Miscellaneous Notes
### LLE Figures
Code to reproduce the LLE panels in Figures 3 and 5 will be provided shortly.
### Figure 7
When viewing the reconstruction of PN trajectories by KCs in PCA space (Figure 7J), be sure to rotate the view to see the reconstruction more clearly. The default view is sometimes not very clear.
### Figure 8
Although data is provided to plot the classification performance of random PN subsets (the blue traces in Figure 8 and S8), the data files are very large and it may take a long time to load them when reproducing these figures. Hence by default, plotting their performance has been disabled. To enable it, set `plotShuffles` to `true` in `code/figures/figure8/MakeFigures.m`. 

# Bugs/Feedback
Please submit any bugs, documentation or feature requests to the [issues tracker](https://github.com/stootoon/ShenTootoonianLaurent2013/issues), or directly by email to <sina.tootoonian@gmail.com>. 

# Downloads

* ShenTootoonianLaurent2013Code [ v0.9](https://dl.dropboxusercontent.com/u/5517909/ShenTootoonianLaurent2013Code_v0.9.zip): Initial public release.
* [Data](https://dl.dropboxusercontent.com/u/5517909/ShenTootoonianLaurent2013Data.zip)