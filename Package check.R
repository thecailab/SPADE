
path <- find.package("SCRIP")
system(paste(shQuote(file.path(R.home("bin"), "R")),
     "CMD", "Rd2pdf", shQuote(path)))

setwd("C:\\Users\\fqin\\Downloads\\SCRIP-main")

devtools::document()
devtools::load_all()
usethis::use_vignette("SCRIPsimu")
usethis::use_testthat()
devtools::test()
usethis::use_news_md()

devtools::build_vignettes()
devtools::check()

devtools::build()

# Useful before submitting to cran
rhub::check_for_cran()


system("E:/DB/Dropbox/Qinfei/Research/SC CNV/Code/FLCNA-package/src/flowCalcCpp.cpp")
library(Rcpp)
compileAttributes()


setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\SC CNV\\Code\\FLCNA-package")
devtools::document()
devtools::build_manual(path=getwd())



setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\Code\\SPADE_package")
devtools::document()
devtools::build_manual(path=getwd())
## add 
# useDynLib(FLCNA, .registration=TRUE)
# exportPattern("^[[:alpha:]]+")
# importFrom(Rcpp, evalCpp)
## in NAMESPACE file





