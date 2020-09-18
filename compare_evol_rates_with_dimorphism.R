## This script was used to analyze the difference in male and female evolutionary rates in
## relation to the level of dichromatism along the same branches.

library(RRphylo)
library(tidyverse)
library(ape)
library(tidytree)

centroids <- data.table::fread('centroids.csv')
tree <- read.nexus('tree.nex')

magnitude <- function(x) {
  if (is.matrix(x)) sqrt(rowSums(x^2)) else sqrt(sum(x^2))
}

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

### Estimate rates using ridge regression
# This can take a long time, so I supply an .rds file instead:
#library(RRphylo)
#RR <- RRphylo(tree, m)
RR <- read_rds('RR.rds')

# Get the estimated ancestral states, combined with tip values
aces_f <- rbind(RR$aces[, 1:3], m[, 1:3])
aces_m <- rbind(RR$aces[, 4:6], m[, 4:6])

# dichromatism:
dE <- sqrt(rowSums((aces_f - aces_m)^2))
# evolutionary rates:
rF <- magnitude(RR$multiple.rates[, 1:3])
rM <- magnitude(RR$multiple.rates[, 4:6])

# Function to obtain the correlation between dichromatism and evolutionary rates.
extract_cors <- function(m, rr, plot = FALSE) {
  aces_f <- rbind(rr$aces[, 1:3], m[, 1:3])
  aces_m <- rbind(rr$aces[, 4:6], m[, 4:6])
  dE <- magnitude(aces_f - aces_m)
  rF <- magnitude(rr$multiple.rates[, 1:3])
  rM <- magnitude(rr$multiple.rates[, 4:6])
  
  return(c(
    ratio  = coef(lm(log(rM / rF) ~ dE))[2],
    ratio_intercept = coef(lm(log(rM / rF) ~ dE))[1]
  ))
}
# This file is generated in compare_evol_rates.R
perms <- read_rds('evolrate_RR_perms.rds')

# observed slope
O_ratio <- extract_cors(m, RR)[3]
# expected slopes
E_ratio <- map_dbl(perms, ~extract_cors(m, .)[3])
# p-value
p_ratio <- pmin(mean(O_ratio >= E_ratio), mean(O_ratio <= E_ratio)) * 2
if (p_ratio == 0) p_ratio <- 0.001

# Make figure 3:
p1d <- data_frame(dE, rF = rF, rM = rM)

p1 <- ggplot(p1d, aes(dE, rM / rF)) +
  geom_point(size = 1, shape = 21, alpha = 0.5, col = 1, fill = 1) + 
  geom_hline(yintercept = 1, color = 'grey60') +
  geom_smooth(method = 'lm', se = FALSE, show.legend = FALSE, col = 'firebrick') +
  scale_y_log10() +
  labs(x = 'dichromatism (||D||)', y = expression(frac('male rate', 'female rate'))) +
  theme_minimal()

xmax <- 63
p1_plus <- p1 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perms, function(.) {
      x <- extract_cors(m, .)
      data.frame(x = 0, y = exp(x[6]), xend = xmax, yend = exp(xmax * x[3] + x[6]))
    } ), 
    size = 0.05, alpha = 0.05
  ) +
  geom_smooth(method = 'lm', se = FALSE, show.legend = FALSE, col = 'firebrick')
ggsave('paper_plots/figure3.png', p1_plus, height = 4, width = 4, dpi = 400)
ggsave('paper_plots/figure3.pdf', p1_plus, height = 4, width = 4)
