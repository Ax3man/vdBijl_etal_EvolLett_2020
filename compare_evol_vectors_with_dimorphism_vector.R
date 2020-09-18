## This script was used to analyze male and female evolutionary vectors, and compare them to the
## direction of dichromatism. This allows us to relate changes in dichromatism to male or female
## specific color evolution.

### Prep -------------------------------------------------------------------------------------------
library(RRphylo)
library(tidyverse)
library(ape)
library(furrr)

centroids <- data.table::fread('centroids.csv')
tree <- read.nexus('tree.nex')

m <- centroids %>% 
  dplyr::select(
    Species, 
    L_mean_females, a_mean_females, b_mean_females,
    L_mean_males, a_mean_males, b_mean_males, 
  ) %>% 
  filter(Species %in% tree$tip.label) %>% 
  as.data.frame() %>% 
  column_to_rownames('Species') %>% 
  as.matrix()
m <- m[tree$tip.label, ]

### Fit and analyze observed -----------------------------------------------------------------------
# Estimate rates using ridge regression
# This can take a long time, so I supply an .rds file instead:
#library(RRphylo)
#RR <- RRphylo(tree, m)
RR <- read_rds('RR.rds')

# set up some helper functions for vector calculations:
magnitude <- function(x) {
  if (is.matrix(x)) sqrt(rowSums(x^2)) else sqrt(sum(x^2))
}
scalar_projection <- function(a, b) { 
  sca_proj <- function(A, B) (A %*% B) / sqrt(sum(B^2))
  if (is.vector(a) & is.vector(b)) return(sca_proj(a, b))
  map_dbl(1:nrow(a), ~sca_proj(a[., ], b[., ]))
}

# Predicted color locations across the tree
aces_f <- rbind(RR$aces[, 1:3], m[, 1:3])
aces_m <- rbind(RR$aces[, 4:6], m[, 4:6])
# The predicted Euclidean distance between male and female phenotypes (i.e. dichromatism)
E <- sqrt(rowSums((aces_m - aces_f)^2))
# Evolutionary vectors in 3d color space
vF <- RR$multiple.rates[, 1:3]
vM <- RR$multiple.rates[, 4:6]
# The vectors from female to male phenotype, i.e. direction and size of dichromatism
vE <- aces_m - aces_f
# The evolutionary rate of dichromatism itself
rE <- sqrt(rowSums(vE ^ 2))
# Evolutionary rates for males and females, i.e. the vector magnitude
rF <- magnitude(RR$multiple.rates[, 1:3])
rM <- magnitude(RR$multiple.rates[, 4:6])
# Effective rates for males and females, i.e. the size of the projection on the the dichromatic direction
rF_eff <- scalar_projection(-vF, vE)
rM_eff <- scalar_projection(vM, vE)
# Vector of change in dichromatism
vdE <- vM - vF
rvdE <- sqrt(rowSums((vdE)^2))
rE_eff <- scalar_projection(vdE, vE)

### Fit and analyze permuted -----------------------------------------------------------------------

perms <- read_rds('evolrate_RR_perms.rds')

extract_cors <- function(rr) {
  aces_f <- rbind(rr$aces[, 1:3], m[, 1:3])
  aces_m <- rbind(rr$aces[, 4:6], m[, 4:6])
  dE <- magnitude(aces_f - aces_m)
  # Estimated evolutionary changes per axis (3d vectors)
  vF <- rr$multiple.rates[, 1:3]
  vM <- rr$multiple.rates[, 4:6]
  # Estimated ancestral state of dichromatism
  vE <- aces_m - aces_f
  # The evolutionary rate of dichromatism
  rE <- sqrt(rowSums(vE ^ 2))
  # Evolutionary rates for males and females, i.e. the vector magnitude
  rF <- magnitude(rr$multiple.rates[, 1:3])
  rM <- magnitude(rr$multiple.rates[, 4:6])
  # Effective rates for males and females, i.e. the size of the projection onto the dichromatic direction
  rF_eff <- scalar_projection(-vF, vE)
  rM_eff <- scalar_projection(vM, vE)
  # Vector of change in dichromatism
  vdE <- vM - vF
  rF_eff2 <- scalar_projection(-vF, vdE)
  rM_eff2 <- scalar_projection(vM, vdE)
  # Effictive rate of change in dichromatism, i.e. change in dichromatism projected onto existing dichromatism
  rE_eff <- scalar_projection(vdE, vE)
  
  m1 <- lm(rM_eff ~ rE_eff)
  m2 <- lm(rF_eff ~ rE_eff)
  
  return(setNames(
    c(coef(m1), coef(m2)),
    c('male_intercept', 'male_slope', 'female_intercept', 'female_slope')
  ))
}

extract_cors(RR)

# observed difference in slopes
O <- diff(extract_cors(RR)[c(2, 4)])
# expected difference in slopes
plan(multisession, workers = 16)
E <- future_map_dbl(perms, ~diff(extract_cors(.)[c(2, 4)]), .progress = TRUE)
plan(sequential)
# p-value
p <- pmin(mean(O >= E), mean(O <= E)) * 2
if (p == 0) p <- 0.001

# Make figure 4
p1 <- ggplot(mapping = aes(rE_eff, rM_eff)) +
  geom_point(alpha = 0.3, shape = 21, fill = 1) +
  theme_minimal()
xmin <- min(rE_eff); xmax <- max(rE_eff)
p1_plus <- p1 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perms, function(.) {
      x <- extract_cors(.)
      data.frame(x = xmin, y = xmin * x[2] + x[1], xend = xmax, yend = xmax * x[2] + x[1])
    } ), 
    size = 0.05, alpha = 0.025
  ) +
  geom_smooth(method = 'lm', se = FALSE, col = 'firebrick') +
  labs(x = 'rate of change in dichromatism', y = 'rate of change in dichromatism\ndue to males')

p2 <- p1 + aes(rE_eff, rF_eff)
p2_plus <- p2 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perms, function(.) {
      x <- extract_cors(.)
      data.frame(x = xmin, y = xmin * x[4] + x[3], xend = xmax, yend = xmax * x[4] + x[3])
    } ), 
    size = 0.05, alpha = 0.025
  ) +
  geom_smooth(method = 'lm', se = FALSE, col = 'firebrick') +
  labs(x = 'rate of change in dichromatism', y = 'rate of change in dichromatism\ndue to females')

cowplot::plot_grid(p1_plus, p2_plus, labels = 'auto')

ggsave('paper_plots/figure4.png', height = 4, width = 7, dpi = 300)
ggsave('paper_plots/figure4.pdf', height = 4, width = 7)
