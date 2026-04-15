# Release Checklist

1. Run `Rscript tools/render_custom_site.R`.
2. Run `Rscript tools/check_release_sync.R`.
3. Run `Rscript -e 'testthat::test_local(stop_on_failure = FALSE)'`.
4. Run `R CMD build .`.
5. Run `R CMD check --no-manual hte3_0.1.1.tar.gz`.
6. Review `NEWS.md`, `README.md`, `docs/`, and `vignettes/` together before tagging.
7. Push the branch and confirm GitHub Actions passes `R-CMD-check`, the custom site sync gate, and the legacy smoke workflow.
8. Merge to `main`, let the site deploy workflow publish `docs/`, and then confirm the GitHub/r-universe release install path.
