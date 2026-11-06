test_that("Test custom term registration", {
  skip_on_cran()
  pkg_path <- file.path(tempdir(), "test_pkg")
  if (dir.exists(pkg_path)) unlink(pkg_path, recursive = TRUE)
  dir.create(pkg_path, recursive = TRUE)

  on.exit({
    if ("iglm.custom" %in% loadedNamespaces()) {
      unloadNamespace("iglm.custom")
    }
    if (dir.exists(pkg_path)) unlink(pkg_path, recursive = TRUE)
    try(remove.packages("iglm.custom", lib = .libPaths()[1]), silent = TRUE)
  }, add = TRUE)

  create_userterms_skeleton(path = pkg_path, pkg_name = "iglm.custom")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(file.path(pkg_path, "iglm.custom"))

  r_cmd <- file.path(R.home("bin"), "R")
  libs_env <- paste(.libPaths(), collapse = .Platform$path.sep)
  status <- system2(
    r_cmd,
    c("CMD", "INSTALL", paste0("--library=", .libPaths()[1]), "."),
    env = c(paste0("R_LIBS=", libs_env), paste0("R_LIBS_USER=", libs_env))
  )
  expect_equal(status, 0)

  library(iglm.custom, lib.loc = .libPaths()[1])
  data(state_twitter)

  res <- statistics(state_twitter[[1]] ~ mutual + my_mutual)
  expect_equal(as.numeric(res[1]), as.numeric(res[2]))
})
