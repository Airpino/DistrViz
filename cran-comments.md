## Submission

This is the first submission of DistrViz.

DistrViz imports DDAWASS (discrete distribution data analysis with the
Wasserstein distance, same maintainer), available on CRAN since [DATE].

## Test environments

* local Windows 11, R 4.5.x
* win-builder: R-devel and R-release
* R-hub: linux (R-devel), macos (R-release), windows (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission.

## Notes for the reviewers

* The interactive versions use the suggested packages 'ggiraph' (with
  'htmlwidgets') and 'plotly'; examples and vignettes run them only if they
  are installed.
* The vignettes use the 'univer' data of the suggested package 'CUB'; without
  it, the code is shown but not run.
