html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

load_release_metadata <- function(root = ".") {
  env <- new.env(parent = baseenv())
  sys.source(file.path(root, "R", "release_metadata.R"), envir = env)
  env
}

format_constructor_call <- function(spec) {
  args <- spec$constructor_args
  if (length(args) == 0L) {
    return(sprintf("%s$new()", spec$learner_class))
  }

  formatted_args <- vapply(
    names(args),
    function(name) paste0(name, " = ", paste(deparse(args[[name]]), collapse = "")),
    character(1)
  )
  sprintf("%s$new(%s)", spec$learner_class, paste(formatted_args, collapse = ", "))
}

render_install_block_html <- function(meta) {
  command <- meta$hte3_install_command()
  lines <- c(
    'if (!requireNamespace("remotes", quietly = TRUE)) {',
    "  install.packages(\"remotes\")",
    "}",
    "",
    command
  )
  sprintf("<pre><code>%s</code></pre>", paste(html_escape(lines), collapse = "\n"))
}

render_default_stack_note_html <- function() {
  paste(
    "<p class=\"panel-note\">",
    "If you do not pass learners manually, <code>hte3</code> uses",
    "<code>get_autoML()</code> for the default nuisance and base-learner stack.",
    "That stack always includes <code>Lrnr_glmnet</code> and <code>Lrnr_gam</code>,",
    "and it adds the <code>earth</code>, <code>ranger</code>, and <code>xgboost</code>",
    "learners when those optional runtime packages are installed.",
    "</p>"
  )
}

render_r_defaults_note_html <- function() {
  paste(
    "<p>",
    "<code>hte_task()</code> defaults the propensity and outcome learners to",
    "<code>get_autoML()</code>, and the fit wrappers default <code>base_learner</code>",
    "to the same stack. The stack always includes its core learners and",
    "adds <code>earth</code>, <code>ranger</code>, and <code>xgboost</code>",
    "components when those optional packages are available.",
    "</p>"
  )
}

render_wrapper_note_html <- function(meta) {
  cate_methods <- paste(toupper(meta$supported_cate_methods()), collapse = ", ")
  crr_methods <- paste(toupper(meta$supported_crr_methods()), collapse = ", ")

  paste(
    "<strong>Wrapper note:</strong>",
    sprintf("<code>fit_cate()</code> supports %s in the high-level API.", cate_methods),
    sprintf("<code>fit_crr()</code> supports %s.", crr_methods),
    "For continuous-treatment CATE, the supported high-level path is the",
    "R-learner under a partially linear <code>A * tau(X)</code> effect model.",
    "Lower-level CRR classes remain available through the advanced interface.",
    "When <code>modifiers</code> are a strict subset of <code>confounders</code>,",
    "prefer DR, EP, or the default two-stage T-learner for the",
    "<code>V</code>-conditional CATE target; the current R-learner instead",
    "targets an overlap-weighted projection and warns at fit time."
  )
}

render_automl_stack_html <- function(meta) {
  specs <- meta$automl_spec_table()
  core_specs <- Filter(function(spec) is.null(spec$required_package), specs)
  optional_specs <- Filter(function(spec) !is.null(spec$required_package), specs)

  lines <- c(
    "Stack$new(",
    paste0("  ", vapply(core_specs, format_constructor_call, character(1)), ","),
    "  # Added when optional runtime packages are installed:",
    paste0("  ", vapply(optional_specs, format_constructor_call, character(1)), ","),
    ")"
  )

  sprintf("<pre><code>%s</code></pre>", paste(html_escape(lines), collapse = "\n"))
}

site_fragment_map <- function(root = ".") {
  meta <- load_release_metadata(root)
  list(
    "docs/r.html" = list(
      install_block = render_install_block_html(meta),
      default_stack_note = render_default_stack_note_html(),
      r_defaults_note = render_r_defaults_note_html()
    ),
    "docs/index.html" = list(
      wrapper_note = render_wrapper_note_html(meta)
    ),
    "docs/sl3.html" = list(
      automl_stack = render_automl_stack_html(meta)
    )
  )
}

replace_marker_block <- function(text, marker, replacement) {
  pattern <- sprintf(
    "([ \t]*)<!-- BEGIN:%s -->([\\s\\S]*?)\\1<!-- END:%s -->",
    marker,
    marker
  )

  if (!grepl(pattern, text, perl = TRUE)) {
    stop(sprintf("Could not find marker block `%s`.", marker), call. = FALSE)
  }

  match_data <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  indent <- match_data[[2]]
  replacement_block <- sprintf(
    "%s<!-- BEGIN:%s -->\n%s\n%s<!-- END:%s -->",
    indent,
    marker,
    replacement,
    indent,
    marker
  )
  gsub(pattern, replacement_block, text, perl = TRUE)
}

render_site_fragments <- function(root = ".", write = FALSE) {
  fragments <- site_fragment_map(root)
  rendered <- list()

  for (path in names(fragments)) {
    full_path <- file.path(root, path)
    text <- paste(readLines(full_path, warn = FALSE), collapse = "\n")

    for (marker in names(fragments[[path]])) {
      text <- replace_marker_block(text, marker, fragments[[path]][[marker]])
    }

    rendered[[path]] <- text

    if (isTRUE(write)) {
      writeLines(text, full_path)
    }
  }

  rendered
}

read_pkgdown_article_slugs <- function(root = ".") {
  lines <- readLines(file.path(root, "_pkgdown.yml"), warn = FALSE)
  in_articles <- FALSE
  slugs <- character()

  for (line in lines) {
    if (grepl("^articles:\\s*$", line)) {
      in_articles <- TRUE
      next
    }

    if (in_articles && grepl("^[A-Za-z_].*:\\s*$", line)) {
      break
    }

    if (in_articles && grepl("^\\s*-\\s*[A-Za-z0-9][A-Za-z0-9-]*\\s*$", line)) {
      slug <- sub("^\\s*-\\s*", "", line)
      if (!(slug %in% c("title", "navbar"))) {
        slugs <- c(slugs, slug)
      }
    }
  }

  unique(slugs)
}
