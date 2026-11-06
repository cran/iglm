test_that("Undirected network with isolated units estimates, simulates, assesses, and predicts cleanly", {
  n <- 6
  # Edges only between node 1 and node 2; nodes 3-6 are isolates
  adj <- matrix(0, n, n)
  adj[1, 2] <- adj[2, 1] <- 1
  nb <- matrix(1, n, n)
  data_obj <- iglm.data(z_network = adj, neighborhood = nb, directed = FALSE)

  sampler <- sampler.iglm(n_burn_in = 10, n_simulation = 5)
  model <- iglm(
    formula = data_obj ~ edges(mode = "local") + degrees,
    sampler = sampler,
    control = control.iglm(max_it = 15, display_progress = FALSE, return_samples = TRUE)
  )

  # Estimation
  expect_no_error(model$estimate())
  expect_true(is.na(model$coef_degrees["degrees 3", 1]))
  expect_true(is.na(model$coef_degrees["degrees 4", 1]))
  expect_false(is.na(model$coef_degrees["degrees 1", 1]))
  expect_false(is.na(model$coef_degrees["degrees 2", 1]))

  # Simulation with degree-0 NA coefficients preserved
  expect_no_error(model$simulate())
  expect_gte(length(model$results$samples), 5)

  # Assessment
  expect_no_error({
    ass <- model$assess(~ degree_distribution + geodesic_distances_distribution, plot = FALSE)
  })

  # Single and multi model assessment plotting
  m2 <- iglm(
    formula = data_obj ~ edges(mode = "local"),
    sampler = sampler,
    control = control.iglm(max_it = 10, display_progress = FALSE)
  )
  m2$estimate()
  ass2 <- m2$assess(~ degree_distribution + geodesic_distances_distribution, plot = FALSE)

  tmp_pdf <- tempfile(fileext = ".pdf")
  pdf(tmp_pdf)
  expect_no_error(model$results$plot(model_assessment = TRUE))
  expect_no_error(model$results$plot(model_assessment = TRUE, ass2))
  dev.off()
  if (file.exists(tmp_pdf)) unlink(tmp_pdf)

  # Predictions
  expect_no_error({
    p_cond <- model$predict(variant = "conditional", type = c("x", "y", "z"))
  })
  expect_no_error({
    p_marg <- model$predict(variant = "marginal", type = c("x", "y", "z"))
  })
  expect_equal(nrow(p_marg$z), n * n)
})

test_that("Directed network with isolates (in=0, out=0) works under estimation, simulation, assessment, and prediction", {
  n <- 5
  # Directed: 1 -> 2 only; node 1 in=0, node 2 out=0; nodes 3, 4, 5 isolates
  adj <- matrix(0, n, n)
  adj[1, 2] <- 1
  nb <- matrix(1, n, n)
  data_obj <- iglm.data(z_network = adj, neighborhood = nb, directed = TRUE)

  sampler <- sampler.iglm(n_burn_in = 10, n_simulation = 5)
  m_dir <- iglm(
    formula = data_obj ~ edges(mode = "local") + degrees,
    sampler = sampler,
    control = control.iglm(max_it = 10, display_progress = FALSE, return_samples = TRUE)
  )

  expect_no_error(m_dir$estimate())
  expect_true(is.na(m_dir$coef_degrees["outdegrees 2", 1]))
  expect_true(is.na(m_dir$coef_degrees["indegrees 1", 1]))
  expect_true(all(is.na(m_dir$coef_degrees[paste0("outdegrees ", 3:5), 1])))
  expect_true(all(is.na(m_dir$coef_degrees[paste0("indegrees ", 3:5), 1])))

  expect_no_error(m_dir$simulate())
  expect_no_error(m_dir$assess(~ degree_distribution + geodesic_distances_distribution, plot = FALSE))

  tmp_pdf <- tempfile(fileext = ".pdf")
  pdf(tmp_pdf)
  expect_no_error(m_dir$results$plot(model_assessment = TRUE))
  dev.off()
  if (file.exists(tmp_pdf)) unlink(tmp_pdf)

  expect_no_error(m_dir$predict(variant = "conditional", type = c("x", "y", "z")))
  expect_no_error(m_dir$predict(variant = "marginal", type = c("x", "y", "z")))
})

test_that("Completely empty network (all isolated units) handles all operations gracefully", {
  n <- 4
  adj <- matrix(0, n, n)
  nb <- matrix(1, n, n)
  data_obj <- iglm.data(z_network = adj, neighborhood = nb, directed = FALSE)

  sampler <- sampler.iglm(n_burn_in = 10, n_simulation = 5)
  m_empty <- iglm(
    formula = data_obj ~ edges(mode = "local") + degrees,
    sampler = sampler,
    control = control.iglm(max_it = 10, display_progress = FALSE, return_samples = TRUE)
  )

  expect_no_error(m_empty$estimate())
  expect_true(all(is.na(m_empty$coef_degrees)))

  expect_no_error(m_empty$simulate())
  expect_no_error(m_empty$assess(~ degree_distribution + geodesic_distances_distribution, plot = FALSE))

  tmp_pdf <- tempfile(fileext = ".pdf")
  pdf(tmp_pdf)
  expect_no_error(m_empty$results$plot(model_assessment = TRUE))
  dev.off()
  if (file.exists(tmp_pdf)) unlink(tmp_pdf)

  expect_no_error(m_empty$predict(variant = "conditional", type = c("x", "y", "z")))
  expect_no_error(m_empty$predict(variant = "marginal", type = c("x", "y", "z")))
})