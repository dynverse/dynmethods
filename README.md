
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build status](https://travis-ci.org/dynverse/dynmethods.svg?branch=master)](https://travis-ci.org/dynverse/dynmethods)

dynmethods
==========

This package contains wrappers for all of the trajectory inference methods included in the [dynverse](https://www.github.com/dynverse/dynverse) review. These wrappers contain code to translate the data structures produced by any of the methods to the common trajectory model (Fig. 1b in [manuscript](https://www.biorxiv.org/content/early/2018/03/05/276907)). The output of each method is transformed into a common trajectory model using [dynwrap](https://www.github.com/dynverse/dynwrap).

Some methods are directly implemented & wrapped inside of R. Other methods, primarily those implemented in python or other languages are wrapped inside a docker container to avoid dependency issues.

Currently implemented are the following wrappers:

| Name                      | Code                                                                                              | Vignette | Containerised                                       |
|:--------------------------|:--------------------------------------------------------------------------------------------------|:---------|:----------------------------------------------------|
| Angle                     | [R/ti\_angle.R\#9](https://github.com/dynverse/dynmethods/blob/master/R/ti_angle.R#9)             |          |                                                     |
| CellRouter                | [containers/cellrouter](https://github.com/dynverse/dynmethods/blob/master/containers/cellrouter) |          | [Yes](https://hub.docker.com/r/dynverse/cellrouter) |
| DPT                       | [R/ti\_dpt.R\#23](https://github.com/dynverse/dynmethods/blob/master/R/ti_dpt.R#23)               |          |                                                     |
| ElPiGraph cyclic          | [containers/elpicycle](https://github.com/dynverse/dynmethods/blob/master/containers/elpicycle)   |          | [Yes](https://hub.docker.com/r/dynverse/elpicycle)  |
| ElPiGraph.R               | [containers/elpigraph](https://github.com/dynverse/dynmethods/blob/master/containers/elpigraph)   |          | [Yes](https://hub.docker.com/r/dynverse/elpigraph)  |
| ElPiGraph linear          | [containers/elpilinear](https://github.com/dynverse/dynmethods/blob/master/containers/elpilinear) |          | [Yes](https://hub.docker.com/r/dynverse/elpilinear) |
| Embeddr                   | [R/ti\_embeddr.R\#12](https://github.com/dynverse/dynmethods/blob/master/R/ti_embeddr.R#12)       |          |                                                     |
| error                     | [R/ti\_error.R\#8](https://github.com/dynverse/dynmethods/blob/master/R/ti_error.R#8)             |          |                                                     |
| Growing Neural Gas        | [R/ti\_gng.R\#12](https://github.com/dynverse/dynmethods/blob/master/R/ti_gng.R#12)               |          |                                                     |
| GPfates                   | [containers/gpfates](https://github.com/dynverse/dynmethods/blob/master/containers/gpfates)       |          | [Yes](https://hub.docker.com/r/dynverse/gpfates)    |
| GrandPrix                 | [containers/grandprix](https://github.com/dynverse/dynmethods/blob/master/containers/grandprix)   |          | [Yes](https://hub.docker.com/r/dynverse/grandprix)  |
| Identity                  | [R/ti\_identity.R\#8](https://github.com/dynverse/dynmethods/blob/master/R/ti_identity.R#8)       |          |                                                     |
| MATCHER                   | [containers/matcher](https://github.com/dynverse/dynmethods/blob/master/containers/matcher)       |          | [Yes](https://hub.docker.com/r/dynverse/matcher)    |
| MERLot                    | [R/ti\_merlot.R\#13](https://github.com/dynverse/dynmethods/blob/master/R/ti_merlot.R#13)         |          |                                                     |
| MFA                       | [R/ti\_mfa.R\#10](https://github.com/dynverse/dynmethods/blob/master/R/ti_mfa.R#10)               |          |                                                     |
| Mpath                     | [R/ti\_mpath.R\#9](https://github.com/dynverse/dynmethods/blob/master/R/ti_mpath.R#9)             |          |                                                     |
| Ouija                     | [R/ti\_ouija.R\#11](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouija.R#11)           |          |                                                     |
| Ouijaflow                 | [containers/ouijaflow](https://github.com/dynverse/dynmethods/blob/master/containers/ouijaflow)   |          | [Yes](https://hub.docker.com/r/dynverse/ouijaflow)  |
| PAGA                      | [containers/paga](https://github.com/dynverse/dynmethods/blob/master/containers/paga)             |          | [Yes](https://hub.docker.com/r/dynverse/paga)       |
| p-Creode                  | [containers/pcreode](https://github.com/dynverse/dynmethods/blob/master/containers/pcreode)       |          | [Yes](https://hub.docker.com/r/dynverse/pcreode)    |
| Periodic Principal Curves | [R/ti\_periodpc.R\#11](https://github.com/dynverse/dynmethods/blob/master/R/ti_periodpc.R#11)     |          |                                                     |
| PhenoPath                 | [R/ti\_phenopath.R\#11](https://github.com/dynverse/dynmethods/blob/master/R/ti_phenopath.R#11)   |          |                                                     |
| projected PAGA            | [containers/praga](https://github.com/dynverse/dynmethods/blob/master/containers/praga)           |          | [Yes](https://hub.docker.com/r/dynverse/praga)      |
| Pseudogp                  | [R/ti\_pseudogp.R\#12](https://github.com/dynverse/dynmethods/blob/master/R/ti_pseudogp.R#12)     |          |                                                     |
| Random                    | [R/ti\_random.R\#8](https://github.com/dynverse/dynmethods/blob/master/R/ti_random.R#8)           |          |                                                     |
| reCAT                     | [R/ti\_recat.R\#15](https://github.com/dynverse/dynmethods/blob/master/R/ti_recat.R#15)           |          |                                                     |
| SCIMITAR                  | [containers/scimitar](https://github.com/dynverse/dynmethods/blob/master/containers/scimitar)     |          | [Yes](https://hub.docker.com/r/dynverse/scimitar)   |
| SCORPIUS                  | [R/ti\_scorpius.R\#43](https://github.com/dynverse/dynmethods/blob/master/R/ti_scorpius.R#43)     |          |                                                     |
| SCOUP                     | [R/ti\_scoup.R\#16](https://github.com/dynverse/dynmethods/blob/master/R/ti_scoup.R#16)           |          |                                                     |
| SCUBA                     | [containers/scuba](https://github.com/dynverse/dynmethods/blob/master/containers/scuba)           |          | [Yes](https://hub.docker.com/r/dynverse/scuba)      |
| Shuffle                   | [R/ti\_shuffle.R\#10](https://github.com/dynverse/dynmethods/blob/master/R/ti_shuffle.R#10)       |          |                                                     |
| Sincell                   | [R/ti\_sincell.R\#15](https://github.com/dynverse/dynmethods/blob/master/R/ti_sincell.R#15)       |          |                                                     |
| SLICE                     | [R/ti\_slice.R\#9](https://github.com/dynverse/dynmethods/blob/master/R/ti_slice.R#9)             |          |                                                     |
| SLICER                    | [R/ti\_slicer.R\#9](https://github.com/dynverse/dynmethods/blob/master/R/ti_slicer.R#9)           |          |                                                     |
| Slingshot                 | [R/ti\_slingshot.R\#14](https://github.com/dynverse/dynmethods/blob/master/R/ti_slingshot.R#14)   |          |                                                     |
| StemID                    | [R/ti\_stemid.R\#7](https://github.com/dynverse/dynmethods/blob/master/R/ti_stemid.R#7)           |          |                                                     |
| Topslam                   | [containers/topslam](https://github.com/dynverse/dynmethods/blob/master/containers/topslam)       |          | [Yes](https://hub.docker.com/r/dynverse/topslam)    |
| TSCAN                     | [R/ti\_tscan.R\#11](https://github.com/dynverse/dynmethods/blob/master/R/ti_tscan.R#11)           |          |                                                     |
| Wanderlust                | [containers/wanderlust](https://github.com/dynverse/dynmethods/blob/master/containers/wanderlust) |          | [Yes](https://hub.docker.com/r/dynverse/wanderlust) |
| Waterfall                 | [R/ti\_waterfall.R\#8](https://github.com/dynverse/dynmethods/blob/master/R/ti_waterfall.R#8)     |          |                                                     |
| Wishbone                  | [containers/wishbone](https://github.com/dynverse/dynmethods/blob/master/containers/wishbone)     |          | [Yes](https://hub.docker.com/r/dynverse/wishbone)   |
