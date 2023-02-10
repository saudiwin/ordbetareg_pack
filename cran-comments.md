## R CMD check results

0 errors | 0 warnings | 1 notes

❯ checking installed package size ... NOTE
    installed size is 13.4Mb
    sub-directories of 1Mb or more:
      data  10.0Mb
      doc    3.3Mb
      
The data sub-directory contains fitted models that are used to illustrate function examples because fitting the model can take up to an hour. 

Submitting to win-release results in 1 ERROR: 

Error(s) in re-building vignettes:
--- re-building 'package_introduction.Rmd' using rmarkdown
Quitting from lines 247-258 (package_introduction.Rmd) 
Error: processing vignette 'package_introduction.Rmd' failed with diagnostics:
could not find function "avg_slopes"
--- failed re-building 'package_introduction.Rmd'

SUMMARY: processing the following file failed:
  'package_introduction.Rmd'

Error: Vignette re-building failed.
Execution halted

Submitting to win-devel results in OK status - no NOTEs or WARNINGs.

The error in win-release appears to be due to the fact the marginaleffects package, which has the avg_slopes function, does not have a binary available for windows for the current version of the package, 0.9.0. This is the required version listed under the Suggests field in ordbetareg but win-release is installing an older version that does not have the avg_slopes function. However, it is necessary to update to 0.9.0 as otherwise the package vignette fails to compile, and I received an email from the CRAN maintainer requesting that I update the package to handle marginaleffects 0.9.0--see build logs at https://cran.r-project.org/web/checks/check_results_ordbetareg.html. If I do not update the package by 2-16-2023, it will be removed from CRAN.

As such, it would seem that this error in win-release should be ignored as it is not clear why it is not installing the current version of marginaleffects, which has been available since Feb 1 2023.

