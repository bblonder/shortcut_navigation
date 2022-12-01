# README

# Author
Benjamin Wong Blonder

# Contact
benjamin.blonder@berkeley.edu

# Last updated
28 November 2022

# Documentation for
Blonder et al., Navigation between initial and desired community states using shortcuts.

# See also
Code used to generate the 'outputs' folder, as well as explanation of file structure, is in: 
https://github.com/michaelhlim/CoexistenceControl.jl

# Instructions
To replicate all analyses in R:
1. source('1 prep datasets.R'), which generates the 'outputs' subfolder.
2. Then source('2 analysis.R'), which creates figures and tables in the 'figures' subfolder.
3. Figure 1 is created separately and appears in the 'figure 1' subfolder.

# Code verified to run in:
R version 4.2.0 (2022-04-22)Platform: x86_64-apple-darwin17.0 (64-bit)Running under: macOS Monterey 12.6.1Matrix products: defaultLAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dyliblocale:[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8attached base packages:[1] stats     graphics  grDevices utils     datasets  methods   base     other attached packages: [1] caret_6.0-93      lattice_0.20-45   forcats_0.5.2     stringr_1.4.1     purrr_0.3.4       readr_2.1.2       tibble_3.1.8      [8] tidyverse_1.3.2   pdp_0.8.1         ranger_0.14.1     colorRamps_2.3.1  igraph_1.3.4      data.table_1.14.2 ggpubr_0.4.0     [15] ggplot2_3.3.6     tidyr_1.2.1       dplyr_1.0.10     loaded via a namespace (and not attached): [1] nlme_3.1-159         fs_1.5.2             lubridate_1.8.0      RColorBrewer_1.1-3   httr_1.4.4           tools_4.2.0          [7] backports_1.4.1      utf8_1.2.2           R6_2.5.1             rpart_4.1.16         DBI_1.1.3            colorspace_2.0-3    [13] nnet_7.3-17          withr_2.5.0          tidyselect_1.1.2     compiler_4.2.0       cli_3.4.0            rvest_1.0.3         [19] xml2_1.3.3           scales_1.2.1         proxy_0.4-27         digest_0.6.29        pkgconfig_2.0.3      parallelly_1.32.1   [25] dbplyr_2.2.1         rlang_1.0.5          readxl_1.4.1         rstudioapi_0.14      generics_0.1.3       jsonlite_1.8.0      [31] ModelMetrics_1.2.2.2 car_3.1-0            googlesheets4_1.0.1  magrittr_2.0.3       Matrix_1.5-1         Rcpp_1.0.9          [37] munsell_0.5.0        fansi_1.0.3          abind_1.4-5          lifecycle_1.0.2      stringi_1.7.8        pROC_1.18.0         [43] carData_3.0-5        MASS_7.3-58.1        plyr_1.8.7           recipes_1.0.1        grid_4.2.0           parallel_4.2.0      [49] listenv_0.8.0        crayon_1.5.1         haven_2.5.1          splines_4.2.0        hms_1.1.2            pillar_1.8.1        [55] ggsignif_0.6.3       future.apply_1.9.1   reshape2_1.4.4       codetools_0.2-18     stats4_4.2.0         reprex_2.0.2        [61] glue_1.6.2           modelr_0.1.9         vctrs_0.4.1          tzdb_0.3.0           foreach_1.5.2        cellranger_1.1.0    [67] gtable_0.3.1         future_1.28.0        assertthat_0.2.1     gower_1.0.0          prodlim_2019.11.13   broom_1.0.1         [73] e1071_1.7-11         rstatix_0.7.0        class_7.3-20         survival_3.4-0       googledrive_2.0.0    gargle_1.2.1        [79] timeDate_4021.104    iterators_1.0.14     hardhat_1.2.0        lava_1.6.10          globals_0.16.1       ellipsis_0.3.2      [85] ipred_0.9-13   