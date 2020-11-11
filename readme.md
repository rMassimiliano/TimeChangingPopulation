Inference in response-adaptive clinical trials when the enrolled population varies over time
============================================================================================

This repository contains an R package including code used to generate
and analyze clinical trials in the paper *Inference in response-adaptive
clinical trials when the enrolled population varies over time*. The manuscript is available upon request.

The package <tt>TrendUtilities</tt> can be installed with the following
code

``` r
R CMD build TrendUtilities
R CMD INSTALL TrendUtilities -l your_path
```

or alternatively using <tt>devtools</tt>.

The file <tt>generata\_data.R</tt> incude steps to generate clinical
trials used for the simulations.
