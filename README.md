# Temporal Clustering by Affinity Propagation (TCAP)

## The Purpose of the Algorithm

A number of different clustering approaches capable of handling time course data exist, with examples including SplineCluster [(Heard et al., 2005)][heard2005] and BHC [(Cooke et al., 2011)][cooke2011]. Typically, such methods will use some form of assessing co-expression, be it through a direct distance metric or comparing the resulting parameters of the gene expression being modelled in some way. However, this leaves open the query of identifying the upstream regulators responsible for the observed co-expression. TCAP sets out to answer this query by using a complex, information-rich distance measure when performing its clustering, which captures a complete regulatory interaction package within its clusters instead of merely co-expression.

## Algorithm Availability outside iPlant

The version of TCAP discussed in [Kiddle et al., 2010][kiddle2010] can be downloaded [here][tcapdownload] under the name TCAP. This is a reimplementation of that algorithm, and due to Matlab being heavily optimised for matrix use the original version is a bit times faster. If your dataset is enormous (~10.000 genes) you may want to consider moving to the Matlab implementation for performance efficiency, but the version provided on iPlant should be able to handle the data just as well if given enough time.

## Basic Input/Output

TCAP accepts a simple expression CSV on input, with gene names in the first column and numerical time points in the first row. Replicates are not supported, average out all your replicates into a single mean value per time point per gene.

The output is a list of genes in each cluster, corresponding expression plots (with the potential to highlight the centroid profile which acts a cluster centre) and input for functional follow-up via BiNGO/MEME.

## How Does it Work?

TCAP combines two components - an information-rich distance measure and a fitting clustering algorithm which accounts for the nature of the metric. The Qian similarity measure [(Qian et al., 2001)][qian2001] allows for the detection of complicated regulatory patterns, such as time shifts and inversions, capturing a richer regulatory interaction within clusters. However, it comes with the drawback of making it impossible for a clustering algorithm to fall back on a simple measure, such as the mean, to evaluate the information within a cluster and use it to make decisions. Affinity Propagation [(Frey and Dueck, 2007)][frey2007] accounts for the complex nature of the distance metric, and extracts clusters out of a similarity matrix in a quick and effective manner. The end result is a partition of the data, with each gene being assigned to a single cluster.

## Test Run

If you want to get a feel for the input and output formatting of TCAP, you can run the program using a minimal demo dataset provided under `cyverseuk/tcap_testdata/yeast_example_10.csv` in Community Data. You can leave all the parameters as defaults.

## Input in Detail

### Expression CSV File

**Mandatory input.** The data that you're going to cluster with TCAP. The first column should contain gene IDs, whilst the first row should contain the time points in numerical form (without any units or any other form of non-numerical information included in the column name). If you have multiple replicates of your time series experiment, average them out to a single value per gene per time point. The data is standardised (scaled to a normal distribution with mean 0 and standard deviation 1) within the algorithm as it's needed for the distance metric to work.

### Self-Similarity Score

**Default:** 0

When computing the Qian similarity measure, some score needs to be assigned to any individual gene compared with itself. Setting this value to 0 will let the algorithm pick out the median of the computed non-identical-profile scores and assign it to every self-comparison to aid the clustering algorithm. A non-zero value will override the median computation with the provided value. It is recommended to leave this parameter at zero.

### Affinity Propagation Iterations

**Default:** 1000

The number of clustering iterations to perform. Altering this number will linearly affect the run time of the portion of the script performing the clustering.

### Affinity Propagation Convergence

**Default:** 100

If this many consecutive iterations of the clustering yield the same final grouping, then the algorithm is deemed to have converged prematurely and the resulting clusters are deemed to be the final output.

### Prior Iteration Weighting

**Default:** 0.9

When computing new affinity propagation matrix values, a safety measure is put in place to avoid numerical oscillations by making the next step's values being a weighted sum of the value at the prior step and the newly computed value. This parameter is the weight given to the previous step's value when computing this sum. Values below 0.5 and above 0.9 are not recommended.

### Highlight Gene List

Given the complex regulatory nature of the TCAP clusters, it may be desired to highlight certain genes (such as transcription factors) in the output for increased ease of result interpretation. Provide a potential gene list with one gene identifier per line, with the gene identifiers matching the ones used in the expression CSV you provided.

### Highlight Centroids in Cluster Plots

Part of TCAP's output is a set of plots, one for each cluster. Seeing how the Qian similarity measure allows for capturing complex regulatory interactions within a single cluster, highlighting the centroid profile (the one acting as the cluster centre to which all other profiles in the cluster pointed to) can help put the temporal shifts, inversions etc. in some form of context. If checked, the centroid will be plotted in red, and its name will be mentioned in the title of the cluster plot. If unchecked, the centroid will not receive any special treatment and it won't be mentioned in the cluster plot title.

## Output in Detail

### `clusters.txt`

The basic cluster export file. Each cluster is denoted by a header line, and then one line per gene in that cluster. Each gene appears in exactly one cluster.

### `plots/`

A folder with .eps-format visualisations of the expression of the genes in each cluster. Since TCAP uses a complex similarity measure, the genes within a cluster won't merely show co-expression, but also other complex regulatory behaviours. A plot is likely to capture a group of genes which form regulatory relationships and not merely co-regulation by some upstream process not identified as part of a cluster.

### `functional_analysis_inputs/`

The resulting TCAP clusters, as exported in basic form in `clusters.txt`, reformatted into formats accepted on input by BiNGO and MEME to make follow-up biological data mining immediately available.

[heard2005]: http://www.pnas.org/content/102/47/16939.short
[cooke2011]: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-399
[qian2001]: http://www.sciencedirect.com/science/article/pii/S0022283600952197?np=y
[kiddle2010]: https://bioinformatics.oxfordjournals.org/content/26/3/355.full
[tcapdownload]: http://www.wsbc.warwick.ac.uk/stevenkiddle/tcap.html
[frey2007]: http://science.sciencemag.org/content/315/5814/972