## R CMD check results

0 errors | 0 warnings | 1 notes

> checking installed package size ... NOTE
    installed size is  5.7Mb
    sub-directories of 1Mb or more:
      data   5.2Mb
      
The data sub-directory contains fitted models that are used to illustrate function examples because fitting the model can take up to an hour. 

Submitting to win-release results in OK status - no NOTEs or WARNINGs.

Submitting to win-devel results in the following WARNING:

Found the following significant warnings:
  Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
See 'd:/RCompile/CRANguest/R-devel/ordbetareg.Rcheck/00install.out' for details.

The ggplot2 3.4.0 package is not currently available in CRAN. To address this issue, I added code that will use the updated function argument if ggplot2 3.4.0 is detected (see plot.R line number 154):

```
if(packageVersion('ggplot2')=="3.4.0") {

  cont_plot <- cont_plot + geom_density(aes(x=true),linewidth=2,colour="gray",alpha=0.7)

} else {

  cont_plot <- cont_plot + geom_density(aes(x=true),size=2,colour="gray",alpha=0.7)

}
```
