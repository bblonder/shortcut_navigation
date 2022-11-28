README

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

## Code verified to run in:
> sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_0.5.2     stringr_1.4.1     purrr_0.3.4       readr_2.1.2       tibble_3.1.8     
 [6] tidyverse_1.3.2   pdp_0.8.1         ranger_0.14.1     colorRamps_2.3.1  igraph_1.3.4     
[11] data.table_1.14.2 ggpubr_0.4.0      ggplot2_3.3.6     tidyr_1.2.1       dplyr_1.0.10     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9          lubridate_1.8.0     lattice_0.20-45     digest_0.6.29       assertthat_0.2.1   
 [6] foreach_1.5.2       utf8_1.2.2          R6_2.5.1            cellranger_1.1.0    backports_1.4.1    
[11] reprex_2.0.2        httr_1.4.4          pillar_1.8.1        rlang_1.0.5         googlesheets4_1.0.1
[16] readxl_1.4.1        rstudioapi_0.14     car_3.1-0           Matrix_1.5-1        labeling_0.4.2     
[21] textshaping_0.3.6   googledrive_2.0.0   munsell_0.5.0       broom_1.0.1         compiler_4.2.0     
[26] modelr_0.1.9        systemfonts_1.0.4   pkgconfig_2.0.3     tidyselect_1.1.2    codetools_0.2-18   
[31] fansi_1.0.3         viridisLite_0.4.1   crayon_1.5.1        tzdb_0.3.0          dbplyr_2.2.1       
[36] withr_2.5.0         grid_4.2.0          jsonlite_1.8.0      gtable_0.3.1        lifecycle_1.0.2    
[41] DBI_1.1.3           magrittr_2.0.3      scales_1.2.1        cli_3.4.0           stringi_1.7.8      
[46] carData_3.0-5       farver_2.1.1        ggsignif_0.6.3      fs_1.5.2            xml2_1.3.3         
[51] ragg_1.2.2          ellipsis_0.3.2      generics_0.1.3      vctrs_0.4.1         cowplot_1.1.1      
[56] RColorBrewer_1.1-3  iterators_1.0.14    tools_4.2.0         glue_1.6.2          hms_1.1.2          
[61] abind_1.4-5         colorspace_2.0-3    gargle_1.2.1        rstatix_0.7.0       rvest_1.0.3        
[66] haven_2.5.1 