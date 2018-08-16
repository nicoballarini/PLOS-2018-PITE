This project contains the article and code for: 

# Subgroup Identification Clinical Trials via the Predicted Individual Treatment Effect


* `/man` Documentation files for functions

* `/paper` Contains R code for figures, tables and supplementary material

* `/R`  R Code for functions in the package

* `/sim` Scenarios considered in the article

Final version of article: <...>


___

First build and install the package in your environment.

Workflow for building the project:
On a terminal, move to the folder and do

    make sim_normal
    make sim_survival
    make figs
    make tables
    make pdf_with_diff
    make supplementary
    make submission

Please note that the simulations may take a long time to compute. To run just one
case, just use nohup `R CMD BATCH --vanilla simulate.R` in e.g. sim/normal/biom10


Session info:

    > sessionInfo(package = NULL)
    R version 3.4.2 (2017-09-28)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 14.04.5 LTS
    
    Matrix products: default
    BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
    LAPACK: /usr/lib/lapack/liblapack.so.3.0
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    [8] base     
    
    other attached packages:
     [1] bindrcpp_0.2             monomvn_1.9-7            lars_1.2                
     [4] pls_2.6-0                ggplot2_2.2.1            PITE_0.1.0              
     [7] R.utils_2.5.0            R.oo_1.22.0              R.methodsS3_1.7.1       
    [10] selectiveInference_1.2.2 survival_2.41-3          intervals_0.15.1        
    [13] knitr_1.20               MASS_7.3-50              glmnet_2.0-13           
    [16] foreach_1.4.3            Matrix_1.2-14            dplyr_0.7.4             
    [19] broom_0.4.2             
    
    loaded via a namespace (and not attached):
     [1] Rcpp_0.12.11     pryr_0.1.2       highr_0.6        compiler_3.4.2  
     [5] plyr_1.8.4       bindr_0.1.1      iterators_1.0.9  tools_3.4.2     
     [9] digest_0.6.12    gtable_0.2.0     tibble_1.3.3     nlme_3.1-131    
    [13] lattice_0.20-35  pkgconfig_2.0.1  rlang_0.1.2      psych_1.6.6     
    [17] mvtnorm_1.0-7    stringr_1.2.0    rprojroot_1.3-2  grid_3.4.2      
    [21] glue_1.2.0       R6_2.2.2         tidyr_0.8.0      reshape2_1.4.2  
    [25] purrr_0.2.4      magrittr_1.5     backports_1.1.2  scales_0.4.1    
    [29] codetools_0.2-15 splines_3.4.2    assertthat_0.2.0 mnormt_1.5-5    
    [33] colorspace_1.3-2 labeling_0.3     quadprog_1.5-5   stringi_1.1.5   
    [37] lazyeval_0.2.1   munsell_0.4.3   





The randomized lasso functions are implemented in the selectiveInference package
in python. We call python from inside R using the `PythonInR` R package. 
To make it work, some steps need to be followed: 
First install python developer files in a terminal.

    sudo apt-get install python-dev

Install PythonInR package in R:

    install.packages("PythonInR")

Install the python selective package from:
https://github.com/selective-inference/Python-software
I followed this on my terminal:

    mkdir script

First install regreg package

    git clone https://github.com/regreg/regreg
    cd regreg
    git checkout b8205eea21fdba690890768a1f181c6b29f0f194
    sudo python setup.py install
    cd ..

Now the selection package:  

    git clone https://github.com/selective-inference/Python-software
    cd Python-software/

Make sure we are using the same version:

    git checkout fc24a63fdc34cffc1f072ce7d6e94f82e5b3402f
  
Now clone the C code that the package uses:

    git clone https://github.com/selective-inference/C-software/
    cd C-software
    git checkout 851279ffb326b145d00af45b87e7d857e3941ec9
    cd ..

Now install the package

    sudo python setup.py install
  
It may be needed to install additional packages in python. 
In my environment I had to install pip 
https://pip.pypa.io/en/stable/installing/#install-pip
and then:

    python -m pip install timeout-decorator --user
    python -m pip install setuptools --user
    python -m pip install cython --user
    python -m pip install mpmath --user
    python -m pip install pandas --user
    python -m pip install nose --user
  



