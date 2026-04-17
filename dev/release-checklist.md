# Release Checklist

1. Run `Rscript -e 'devtools::document()'` if roxygen comments changed.
2. Run `Rscript -e 'pkgdown::build_site()'`.
3. Run `Rscript -e 'testthat::test_local(stop_on_failure = FALSE)'`.
4. Run `R CMD build .`.
5. Run `R CMD check --no-manual hte3_0.1.1.tar.gz`.
6. Review `NEWS.md`, `README.md`, `vignettes/`, and generated `docs/` output together before tagging.
7. Push the branch and confirm GitHub Actions passes `R-CMD-check`, the site deploy workflow, and the legacy smoke workflow.
8. Merge to `main`, let the site deploy workflow publish `docs/`, and then confirm the GitHub/r-universe release install path.
