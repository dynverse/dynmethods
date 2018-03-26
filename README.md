
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build status](https://travis-ci.org/dynverse/dynmethods.svg?branch=master)](https://travis-ci.org/dynverse/dynmethods)

dynmethods
==========

This package contains wrappers for all of the trajectory inference methods included in the [dynverse](https://www.github.com/dynverse/dynverse) review. These wrappers contain code to translate the data structures produced by any of the methods to the common trajectory model (Fig. 1b in [manuscript](https://www.biorxiv.org/content/early/2018/03/05/276907)). There are several common post-processing functions provided by the [dynwrap](https://www.github.com/dynverse/dynwrap) package (Supp. Fig. 19b in [manuscript](https://www.biorxiv.org/content/early/2018/03/05/276907)), which the different wrappers are allowed to make use of.

Currently implemented are the following wrappers:

| Name                 | Code                                                                                            | Vignette                                                                         |
|:---------------------|:------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------|
| AGA                  | [R/ti\_aga.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_aga.R#L3)             | [aga.md](https://github.com/dynverse/dynmethods/blob/master/inst/doc/aga.md)     |
| AGA pseudotime       | [R/ti\_aga.R\#L7](https://github.com/dynverse/dynmethods/blob/master/R/ti_aga.R#L7)             | [agapt.md](https://github.com/dynverse/dynmethods/blob/master/inst/doc/agapt.md) |
| Arc-tangent          | [R/ti\_atan.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_atan.R#L3)           |                                                                                  |
| cellTree with gibbs  | [R/ti\_celltree.R\#L7](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L7)   |                                                                                  |
| cellTree with maptpx | [R/ti\_celltree.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L3)   |                                                                                  |
| cellTree with vem    | [R/ti\_celltree.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L11) |                                                                                  |
| Component 1          | [R/ti\_comp1.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_comp1.R#L3)         |                                                                                  |
| Control: identity    | [R/ti\_identity.R\#L6](https://github.com/dynverse/dynmethods/blob/master/R/ti_identity.R#L6)   |                                                                                  |
| Control: manual      | [R/ti\_manual.R\#L6](https://github.com/dynverse/dynmethods/blob/master/R/ti_manual.R#L6)       |                                                                                  |
| Control: random      | [R/ti\_random.R\#L6](https://github.com/dynverse/dynmethods/blob/master/R/ti_random.R#L6)       |                                                                                  |
| Control: shuffle     | [R/ti\_shuffle.R\#L6](https://github.com/dynverse/dynmethods/blob/master/R/ti_shuffle.R#L6)     |                                                                                  |
| DPT                  | [R/ti\_dpt.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_dpt.R#L3)             |                                                                                  |
| Embeddr              | [R/ti\_embeddr.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_embeddr.R#L3)     |                                                                                  |
| GPfates              | [R/ti\_gpfates.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_gpfates.R#L3)     |                                                                                  |
| Growing Neural Gas   | [R/ti\_gng.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_gng.R#L3)             |                                                                                  |
| mfa                  | [R/ti\_mfa.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_mfa.R#L3)             |                                                                                  |
| Monocle DDRTree      | [R/ti\_monocle.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_monocle.R#L3)     |                                                                                  |
| Monocle ICA          | [R/ti\_monocle.R\#L7](https://github.com/dynverse/dynmethods/blob/master/R/ti_monocle.R#L7)     |                                                                                  |
| Mpath                | [R/ti\_mpath.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_mpath.R#L3)         |                                                                                  |
| ouija                | [R/ti\_ouija.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouija.R#L3)         |                                                                                  |
| ouijaflow            | [R/ti\_ouijaflow.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouijaflow.R#L3) |                                                                                  |
| Periodic PrinCurve   | [R/ti\_periodpc.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_periodpc.R#L3)   |                                                                                  |
| PhenoPath            | [R/ti\_phenopath.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_phenopath.R#L3) |                                                                                  |
| pseudogp             | [R/ti\_pseudogp.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_pseudogp.R#L3)   |                                                                                  |
| recat                | [R/ti\_recat.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_recat.R#L3)         |                                                                                  |
| SCIMITAR             | [R/ti\_scimitar.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_scimitar.R#L3)   |                                                                                  |
| SCORPIUS             | [R/ti\_scorpius.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_scorpius.R#L3)   |                                                                                  |
| SCOUP                | [R/ti\_scoup.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_scoup.R#L3)         |                                                                                  |
| SCUBA                | [R/ti\_scuba.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_scuba.R#L3)         |                                                                                  |
| Sincell              | [R/ti\_sincell.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_sincell.R#L3)     |                                                                                  |
| SLICE                | [R/ti\_slice.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_slice.R#L3)         |                                                                                  |
| SLICER               | [R/ti\_slicer.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_slicer.R#L3)       |                                                                                  |
| Slingshot            | [R/ti\_slingshot.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_slingshot.R#L3) |                                                                                  |
| StemID               | [R/ti\_stemid.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_stemid.R#L3)       |                                                                                  |
| topslam              | [R/ti\_topslam.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_topslam.R#L3)     |                                                                                  |
| TSCAN                | [R/ti\_tscan.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_tscan.R#L3)         |                                                                                  |
| Waterfall            | [R/ti\_waterfall.R\#L3](https://github.com/dynverse/dynmethods/blob/master/R/ti_waterfall.R#L3) |                                                                                  |
