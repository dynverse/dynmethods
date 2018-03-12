
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build status](https://travis-ci.org/dynverse/dynmethods.svg?branch=master)](https://travis-ci.org/dynverse/dynmethods)

dynmethods
==========

This package contains wrappers for all of the trajectory inference methods included in the [dynverse](https://www.github.com/dynverse/dynverse) review. These wrappers contain code to translate the data structures produced by any of the methods to the common trajectory model (Fig 1b in manuscript). There are several common post-processing functions provided by the [dynwrap](https://www.github.com/dynverse/dynwrap) package (Supp. Fig. 19b in manuscript), which the different wrappers are allowed to make use of.

Currently implemented are the following wrappers:

-   [cellTree with maptpx](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L3)
-   [cellTree with gibbs](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L7)
-   [cellTree with vem](https://github.com/dynverse/dynmethods/blob/master/R/ti_celltree.R#L11)
-   [Component 1](https://github.com/dynverse/dynmethods/blob/master/R/ti_compone.R#L3)
-   [DPT](https://github.com/dynverse/dynmethods/blob/master/R/ti_dpt.R#L3)
-   [embeddr](https://github.com/dynverse/dynmethods/blob/master/R/ti_embeddr.R#L3)
-   [GPfates](https://github.com/dynverse/dynmethods/blob/master/R/ti_gpfates.R#L3)
-   [identity](https://github.com/dynverse/dynmethods/blob/master/R/ti_identity.R#L6)
-   [manual](https://github.com/dynverse/dynmethods/blob/master/R/ti_manual.R#L6)
-   [mfa](https://github.com/dynverse/dynmethods/blob/master/R/ti_mfa.R#L3)
-   [monocle with DDRTree](https://github.com/dynverse/dynmethods/blob/master/R/ti_monocle.R#L3)
-   [monocle with ICA](https://github.com/dynverse/dynmethods/blob/master/R/ti_monocle.R#L7)
-   [Mpath](https://github.com/dynverse/dynmethods/blob/master/R/ti_mpath.R#L3)
-   [ouija](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouija.R#L3)
-   [ouijaflow](https://github.com/dynverse/dynmethods/blob/master/R/ti_ouijaflow.R#L3)
-   [phenopath](https://github.com/dynverse/dynmethods/blob/master/R/ti_phenopath.R#L3)
-   [pseudogp](https://github.com/dynverse/dynmethods/blob/master/R/ti_pseudogp.R#L3)
-   [random](https://github.com/dynverse/dynmethods/blob/master/R/ti_random.R#L6)
-   [recat](https://github.com/dynverse/dynmethods/blob/master/R/ti_recat.R#L3)
-   [SCIMITAR](https://github.com/dynverse/dynmethods/blob/master/R/ti_scimitar.R#L3)
-   [SCORPIUS](https://github.com/dynverse/dynmethods/blob/master/R/ti_scorpius.R#L3)
-   [SCOUP](https://github.com/dynverse/dynmethods/blob/master/R/ti_scoup.R#L3)
-   [SCUBA](https://github.com/dynverse/dynmethods/blob/master/R/ti_scuba.R#L3)
-   [shuffle](https://github.com/dynverse/dynmethods/blob/master/R/ti_shuffle.R#L6)
-   [sincell](https://github.com/dynverse/dynmethods/blob/master/R/ti_sincell.R#L3)
-   [SLICE](https://github.com/dynverse/dynmethods/blob/master/R/ti_slice.R#L3)
-   [SLICER](https://github.com/dynverse/dynmethods/blob/master/R/ti_slicer.R#L3)
-   [slingshot](https://github.com/dynverse/dynmethods/blob/master/R/ti_slingshot.R#L3)
-   [StemID](https://github.com/dynverse/dynmethods/blob/master/R/ti_stemid.R#L3)
-   [topslam](https://github.com/dynverse/dynmethods/blob/master/R/ti_topslam.R#L3)
-   [TSCAN](https://github.com/dynverse/dynmethods/blob/master/R/ti_tscan.R#L3)
-   [Waterfall](https://github.com/dynverse/dynmethods/blob/master/R/ti_waterfall.R#L3)
