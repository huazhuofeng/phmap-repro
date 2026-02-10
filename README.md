# phmap-repro
Reproductive analysis code for PH-Map atlas

### Statistical analysis
All statistical computations were executed in R (v4.4.1). Data distribution and variance homogeneity were evaluated using the Shapiro-Wilk (n ≤ 5000) or Anderson-Darling (n > 5000) test and Bartlett’s test, respectively. Parametric data were analyzed via unpaired Student’s t-test (two groups) or one-way/two-way ANOVA followed by Tukey’s HSD post hoc test; Welch’s ANOVA with Games-Howell test was substituted if variances were unequal. Non-parametric datasets were assessed using Mann-Whitney U test or Kruskal-Wallis test with Dunn’s post hoc tests. Data are expressed as mean ± SD unless specified. Statistical significance was set at P < 0.05, with all analyses performed in a blinded manner.

> Code for pipeline can be found at https://github.com/Doctorluka/Multiple_Group_Statistics.

### PH-Map tools
PH-Map: Multi-task cell type classification package for single-cell RNA sequencing data.

> Reproductive code for PH-Map tools can be found at https://github.com/Doctorluka/PH-Map.


### Genome alignment and quantification
Raw reads were uniformly reprocessed using a STARsolo-based pipeline with a consistent reference genome.

> Code for pipeline can be found at https://github.com/jarxunlai/scripts_for_publication.

### Integration benchmarking
To evaluate batch correction and data harmonization, multiple integration algorithms were applied, categorized by their output types: embeddings, full features, or kNN graphs. Integration was primarily conducted in Scanpy, with select methods implemented in Seurat (4.4.1) for compatibility.

> Code for pipeline can be found at https://github.com/Doctorluka/ScibPlot.
