source("tools/site_support.R")

fail <- function(message) {
  stop(message, call. = FALSE)
}

assert_contains <- function(path, pattern, fixed = TRUE) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  found <- if (isTRUE(fixed)) {
    grepl(pattern, text, fixed = TRUE)
  } else {
    grepl(pattern, text, perl = TRUE)
  }

  if (!isTRUE(found)) {
    fail(sprintf("Expected `%s` to contain `%s`.", path, pattern))
  }
}

check_rendered_site_fragments <- function(root = ".") {
  rendered <- render_site_fragments(root = root, write = FALSE)

  for (path in names(rendered)) {
    current <- paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
    if (!identical(rendered[[path]], current)) {
      fail(sprintf("Custom site fragments are out of sync in `%s`. Run `Rscript tools/render_custom_site.R`.", path))
    }
  }
}

check_vignette_inventory <- function(root = ".") {
  slugs <- read_pkgdown_article_slugs(root)
  readme <- file.path(root, "README.md")
  site_page <- file.path(root, "docs", "r.html")

  for (slug in slugs) {
    assert_contains(readme, sprintf("vignettes/%s.Rmd", slug))
    assert_contains(site_page, sprintf("vignettes/%s.Rmd", slug))
  }
}

check_release_docs <- function(root = ".") {
  meta <- load_release_metadata(root)
  install_command <- meta$hte3_install_command()

  assert_contains(file.path(root, "README.md"), install_command)
  assert_contains(file.path(root, "docs", "r.html"), install_command)
  assert_contains(file.path(root, "docs", "site.js"), sprintf('defaultMethod: "%s"', meta$default_cate_method()))
  assert_contains(file.path(root, "docs", "site.js"), sprintf('defaultMethod: "%s"', meta$default_crr_method()))
  assert_contains(file.path(root, "man", "get_autoML.Rd"), "available safe subset")
  assert_contains(file.path(root, "man", "Lrnr_cate_T.Rd"), "first-stage contrast")
  assert_contains(file.path(root, "man", "Lrnr_crr_T.Rd"), "first-stage")
}

check_rendered_site_fragments()
check_vignette_inventory()
check_release_docs()
