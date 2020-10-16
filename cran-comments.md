## Test environments
* local Linux Ubuntu 18.04, R 4.0.3
* win-builder (R-release 4.0.3) via devtools::check_win_release() on 16-Oct-2020
* win-builder (R-devel 2020-10-15 r79342) via devtools::check_win_devel() on 16-Oct-2020
* rhub::check_for_cran() on 16-Oct-2020

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE: 
New submission. 
MixSIAR version 3.1.11 was previously on CRAN and archived on 2020-09-25 as check problems were not corrected in time.
There was an ERROR in MKL: https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/MixSIAR.out
This check error was due to a failed test, which we have corrected.
Status was 'OK' on all 12 systems: https://cran-archive.r-project.org/web/checks/2020/2020-09-25_check_results_MixSIAR.html
I have also updated my maintainer email: 'Brian Stock <bstock09@gmail.com>'

## Downstream dependencies
There were no downstream dependencies for this package as of version 3.1.11.