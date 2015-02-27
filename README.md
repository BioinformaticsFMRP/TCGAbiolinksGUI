# README #

### What is this repository for? ###

* Quick summary: biomics package
* Version: 0.1

### How do I get set up? ###

* Summary of set up
Download the package [here](https://bitbucket.org/biomics/biomics/downloads) and install with the following command:

```r
install.packages(path_to_file, repos = NULL, type="source")
```
* Configuration
* Dependencies
```
R (>= 3.0),
Pacakges:
  shinydashboard (>= 0.2.2)
  shiny (>= 0.11.1),
  downloader,
  gdata,
  RCurl,
  R.utils,
  rjson,
  reutils,
  XML
```
* Database configuration

* Creating documentation
```r
roxygen2::roxygenise() # create man pages from code
```

* How to run tests
```r
devtools::test() # execute unit test
```
* Deployment instructions

### Contribution guidelines ###

* Writing tests
[Guide](http://r-pkgs.had.co.nz/tests.html)
[Paper](http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf)
* Documentation
[Guide](https://support.rstudio.com/hc/en-us/articles/200532317-Writing-Package-Documentation)
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin: Tiago C. Silva (tiagochst [at] gmail.com)
* Other community or team contact: Tiago Mendes (tmendesilva [at] gmail.com)
