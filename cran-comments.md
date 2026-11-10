## Submission

This is the first submission of DistrViz.

DistrViz imports DDAWASS (>= 1.0.0) and HistDAWass (>= 1.0.10), by the same
maintainer; DDAWASS was published on CRAN on 2026-11-04.

## Test environments

* local Windows 11, R 4.5.3
* win-builder: R-devel and R-release
* GitHub Actions: windows-latest, macos-latest (R release), ubuntu-latest (R release and devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission.

## Notes for the reviewers

* The interactive versions use the suggested packages 'ggiraph' (with
  'htmlwidgets') and 'plotly'; the maps in the vignettes use 'maps'.
  Examples and vignettes run them only if they are installed.
* The vignettes use the 'univer' data of the suggested package 'CUB'.
