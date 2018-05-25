
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build status](https://travis-ci.org/dynverse/dynmethods.svg?branch=master)](https://travis-ci.org/dynverse/dynmethods)

dynmethods
==========

This package contains wrappers for all of the trajectory inference methods included in the [dynverse](https://www.github.com/dynverse/dynverse) review. These wrappers contain code to translate the data structures produced by any of the methods to the common trajectory model (Fig. 1b in [manuscript](https://www.biorxiv.org/content/early/2018/03/05/276907)). There are several common post-processing functions provided by the [dynwrap](https://www.github.com/dynverse/dynwrap) package (Supp. Fig. 19b in [manuscript](https://www.biorxiv.org/content/early/2018/03/05/276907)), which the different wrappers are allowed to make use of.

Currently implemented are the following wrappers:

| Name               | Code                                                                                              | Vignette |
|:-------------------|:--------------------------------------------------------------------------------------------------|:---------|
| Angle              | [R/ti\_angle.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_angle.R#L11)         |          |
| Control: error     | [R/ti\_error.R\#L10](https://github.com/dynverse/dynmethods/blob/master/R/ti_error.R#L10)         |          |
| Control: identity  | [R/ti\_identity.R\#L10](https://github.com/dynverse/dynmethods/blob/master/R/ti_identity.R#L10)   |          |
| Control: random    | [R/ti\_random.R\#L10](https://github.com/dynverse/dynmethods/blob/master/R/ti_random.R#L10)       |          |
| Control: shuffle   | [R/ti\_shuffle.R\#L12](https://github.com/dynverse/dynmethods/blob/master/R/ti_shuffle.R#L12)     |          |
| DPT                | [R/ti\_dpt.R\#L25](https://github.com/dynverse/dynmethods/blob/master/R/ti_dpt.R#L25)             |          |
| Embeddr            | [R/ti\_embeddr.R\#L14](https://github.com/dynverse/dynmethods/blob/master/R/ti_embeddr.R#L14)     |          |
| GPfates            | [R/ti\_gpfates.R\#L14](https://github.com/dynverse/dynmethods/blob/master/R/ti_gpfates.R#L14)     |          |
| Growing Neural Gas | [R/ti\_gng.R\#L13](https://github.com/dynverse/dynmethods/blob/master/R/ti_gng.R#L13)             |          |
| MATCHER            | [R/ti\_matcher.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_matcher.R#L11)     |          |
| MERLoT             | [R/ti\_merlot.R\#L15](https://github.com/dynverse/dynmethods/blob/master/R/ti_merlot.R#L15)       |          |
| mfa                | [R/ti\_mfa.R\#L12](https://github.com/dynverse/dynmethods/blob/master/R/ti_mfa.R#L12)             |          |
| Mpath              | [R/ti\_mpath.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_mpath.R#L11)         |          |
| ouija              | [R/ti\_ouija.R\#L13](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouija.R#L13)         |          |
| ouijaflow          | [R/ti\_ouijaflow.R\#L10](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouijaflow.R#L10) |          |
| PAGA               | [R/ti\_paga.R\#L43](https://github.com/dynverse/dynmethods/blob/master/R/ti_paga.R#L43)           |          |
| Periodic PrinCurve | [R/ti\_periodpc.R\#L13](https://github.com/dynverse/dynmethods/blob/master/R/ti_periodpc.R#L13)   |          |
| PhenoPath          | [R/ti\_phenopath.R\#L13](https://github.com/dynverse/dynmethods/blob/master/R/ti_phenopath.R#L13) |          |
| pseudogp           | [R/ti\_pseudogp.R\#L14](https://github.com/dynverse/dynmethods/blob/master/R/ti_pseudogp.R#L14)   |          |
| reCAT              | [R/ti\_recat.R\#L17](https://github.com/dynverse/dynmethods/blob/master/R/ti_recat.R#L17)         |          |
| SCIMITAR           | [R/ti\_scimitar.R\#L15](https://github.com/dynverse/dynmethods/blob/master/R/ti_scimitar.R#L15)   |          |
| SCORPIUS           | [R/ti\_scorpius.R\#L45](https://github.com/dynverse/dynmethods/blob/master/R/ti_scorpius.R#L45)   |          |
| SCOUP              | [R/ti\_scoup.R\#L18](https://github.com/dynverse/dynmethods/blob/master/R/ti_scoup.R#L18)         |          |
| SCUBA              | [R/ti\_scuba.R\#L15](https://github.com/dynverse/dynmethods/blob/master/R/ti_scuba.R#L15)         |          |
| Sincell            | [R/ti\_sincell.R\#L17](https://github.com/dynverse/dynmethods/blob/master/R/ti_sincell.R#L17)     |          |
| SLICE              | [R/ti\_slice.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_slice.R#L11)         |          |
| SLICER             | [R/ti\_slicer.R\#L11](https://github.com/dynverse/dynmethods/blob/master/R/ti_slicer.R#L11)       |          |
| Slingshot          | [R/ti\_slingshot.R\#L15](https://github.com/dynverse/dynmethods/blob/master/R/ti_slingshot.R#L15) |          |
| StemID             | [R/ti\_stemid.R\#L9](https://github.com/dynverse/dynmethods/blob/master/R/ti_stemid.R#L9)         |          |
| topslam            | [R/ti\_topslam.R\#L12](https://github.com/dynverse/dynmethods/blob/master/R/ti_topslam.R#L12)     |          |
| TSCAN              | [R/ti\_tscan.R\#L13](https://github.com/dynverse/dynmethods/blob/master/R/ti_tscan.R#L13)         |          |
| Wanderlust         | [R/ti\_wishbone.R\#L49](https://github.com/dynverse/dynmethods/blob/master/R/ti_wishbone.R#L49)   |          |
| Waterfall          | [R/ti\_waterfall.R\#L10](https://github.com/dynverse/dynmethods/blob/master/R/ti_waterfall.R#L10) |          |
| Wishbone           | [R/ti\_wishbone.R\#L45](https://github.com/dynverse/dynmethods/blob/master/R/ti_wishbone.R#L45)   |          |
