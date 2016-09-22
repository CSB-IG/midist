# chrom_mi_dist - Distribution plot of MI from adjacency matrices.

This script takes an adjacency matrix of Mutual Information (MI) between genes in a pool
encoded as [SIF files](http://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format )
and builds a lineplot showing the probability distribution of MI.

It will search which genes belong to each chromosome,
and plot the distribution of MI values in the following subsets:

-  All the genes in the matrix.
-  Gene pairs belonging to the same chromosome.
-  Gene pairs belonging to different chromosomes.

# Prerequisites

- R
- mk (usually found in `9base` or `plan9port` package)

# Quickstart

1.  Make a directory named `data`:

```
mk init
```

2.  Open `R` in this directory.

3.  If not done automatically, run on your R console:

```
packrat::restore()
q("no")
```

4.  Add SIF files to the data file.

5.  Run the analysis.

```
mk
```

    Your results will be on the `results` folder.
