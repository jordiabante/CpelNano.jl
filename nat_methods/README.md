# Data

## CpelNano Output Format

The suffix of the file indicates the type of information:

- ex: expected value of X ∈ {-1,1}, ie E[X]
- exx: expected value of XX ∈ {-1,1}, ie E[XX]
- mml: mean methylation level μ(X)
- tmml: mean methylation level test statistic
- nme: normalized methylation entropy h(X)
- tnme: normalized methylation entropy test statistic
- theta: parameter vector θ
- tcmd: coefficient of methylation divergence

The columns have the following meaning when the suffix is not tmml, tnme or tcmd:

    | chr | st | end | quantity |

and the following one otherwise

    | chr | st | end | test statistic | p-value | q-value (BH) |

## Content of subfolders

The *real_data* folder contains all the necessary data to replicate the figures and tables in the paper.

── real_data
│   ├── breast_cancer_cas9_data: data used in the targeted analysis performed in the paper.
│   │   ├── cpelnano: output of cpelnano for cas-9 nanopore breast cancer data used in the paper.
│   │   │             Each *cpelnano_* file has the cell line name (MCF10a/MDAMB231), replicate ID (1-10),
│   │   │             and the quantity type (ex, exx, mml, nme, theta).  Each *grp_comp_* file is part of the
│   │   │             differential analysis output for control vs control or control vs cancer as indicated
│   │   │             by the cell lines.
│   │   └── target_regions: contains target regions.
│   └── gm12878_data: output of cpelnano in chr22 of GM12878 data.
│       ├── cpelnano_nanopore: output for nanopore sequencing data and nanopore read length analysis.
│       ├── cpelnano_wgbs: output for WGBS data.

The *simulations* contains all the results obtained in our simulations.

── simulations
    ├── nanopolish_benchmark: benchmarking results for nanopolish. First column contains ground-truth,
    |   second column contains the LRT, and the third column contains the actual call made by nanopolish.
    ├── groundtruth: groundtruth quantities obtained for chromosome 22 from GM12878 WGBS data.
    └── state_of_the_art_comparison: output for CpelNano modeling noise (cpelnano), CpelNano without modeling noise
    (cpelnano_nonoise) and empirical estimates for different levels of noise (sigma=2.0-3.5). 
