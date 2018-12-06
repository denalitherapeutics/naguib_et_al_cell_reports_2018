options(repos = c(CRAN = "http://cran.rstudio.com"),
        Ncpus = parallel::detectCores())
switchr::switchTo("supt4h1")
install.packages("remotes")
remotes::install_github("denalitherapeutics/naguib_et_al_cell_reports_2019",
                        dependencies = FALSE, upgrade = "never")
