parasite - A package to calculate standard parasitology measures.

I created this small R package to streamline basic parasitology calculations. Many people use quantitative parasitology online and it works just fine that way; however, I wanted to be able to code all my analyses for my thesis. Subsequently, when I finished writing my code, many of the analyses that I did were easy to implement into functions and I made a package for it. I do not have plans to publish this package to the CRAN any time soon as I believe it may be too niche to be applicable. However, if given enough interest I will publish it and write a methods paper for it.

To install the package currently, make sure you have the R package 'devtools' installed and run this code verbatim: 

devtools::install_github("coltonp2022/parasite").

This will install the version currently on github. Each function should be documented relatively well and there is an example dataset that can be used. Additionally, all code can be used within the piping system of the tidyverse. 
