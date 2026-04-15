if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Larsvanderlaan/hte3",
  ref = "legacy-paper-repro",
  upgrade = "never",
  dependencies = TRUE
)
