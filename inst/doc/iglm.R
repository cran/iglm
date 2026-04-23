## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  fig.width = 7, 
  fig.height = 5
)
library(iglm)

## -----------------------------------------------------------------------------
n_actor <- 100

attribute_info <- rnorm(n_actor)
attribute_cov <- diag(attribute_info)
edge_cov <- outer(attribute_info, attribute_info, FUN = function(x, y) {
  abs(x - y)
})
set.seed(123)

alpha <- 0.3
block <- matrix(nrow = 50, ncol = 50, data = 1)
neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 50, block, simplify = FALSE)))

overlapping_degree <- 0.5
neighborhood <- matrix(nrow = n_actor, ncol = n_actor, data = 0)
block <- matrix(nrow = 5, ncol = 5, data = 0)
size_neighborhood <- 5
size_overlap <- ceiling(size_neighborhood * overlapping_degree)

end <- floor((n_actor - size_neighborhood) / size_overlap)
for (i in 0:end) {
  neighborhood[(1 + size_overlap * i):(size_neighborhood + size_overlap * i), (1 + size_overlap * i):(size_neighborhood + size_overlap * i)] <- 1
}
neighborhood[(n_actor - size_neighborhood + 1):(n_actor), (n_actor - size_neighborhood + 1):(n_actor)] <- 1

type_x <- "binomial"
type_y <- "binomial"
object <- iglm.data(neighborhood = neighborhood, directed = F, type_x = type_x, type_y = type_y, n_actor = n_actor)

## -----------------------------------------------------------------------------
formula <- object ~ edges + attribute_y + attribute_x + degrees

## -----------------------------------------------------------------------------
# Parameters of edges(mode = "local"), attribute_y, and attribute_x
gt_coef <- c(3, -1, -1)
# Parameters for degree effect
gt_coef_degrees <- c(rnorm(n = n_actor, -2, 1))
# Define the sampler
sampler_tmp <- sampler.iglm(
  n_burn_in = 100, n_simulation = 10,
  sampler_x = sampler.net.attr(n_proposals = n_actor * 10),
  sampler_y = sampler.net.attr(n_proposals = n_actor * 10),
  sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 10),
  init_empty = F
)

model_tmp_new <- iglm(
  formula = formula,
  coef = gt_coef, coef_degrees = gt_coef_degrees, sampler = sampler_tmp,
  control = control.iglm(accelerated = F, max_it = 200, display_progress = F)
)

## -----------------------------------------------------------------------------
# Simulate new networks
model_tmp_new$simulate()
# Get the samples
tmp <- model_tmp_new$get_samples()

## -----------------------------------------------------------------------------
# First set the first simulated network as the target for estimation
model_tmp_new$set_target(tmp[[1]])
model_tmp_new$estimate()
model_tmp_new$iglm.data$degree_distribution(plot = TRUE)

## -----------------------------------------------------------------------------
model_tmp_new$assess(formula = ~ degree_distribution +
  geodesic_distances_distribution + edgewise_shared_partner_distribution + mcmc_diagnostics)
model_tmp_new$results$plot(model_assessment = T)

