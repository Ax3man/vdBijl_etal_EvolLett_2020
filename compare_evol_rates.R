## This script was used to compare the average evolutionary rates of males and females.
## It also forms the basis of some of the other analyses.

### Load packages and data
library(tidyverse)
library(ape)
library(tidytree)
library(ggtree)

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

### Permutate sex labels to get a null distribution:
fit_permuted_RR <- function(m) {
  require(RRphylo)
  # randomly choose whether to keep or swap the male/female label
  r <- sample(c(TRUE, FALSE), nrow(m), replace = TRUE)
  m_swapped <- m[, c(4:6, 1:3)]
  colnames(m_swapped) <- colnames(m)
  
  # create a permuted matrix and reorder
  m_perm <- rbind(m[!r, ], m_swapped[r, ])
  m_perm <- m_perm[tree$tip.label, ]
  RRphylo(tree, m_perm, clus = 0)
}
# This will be slow. I have the .rds file for perms, but it's too large for github.
# This code will reproduce it, but you can also shoot me an email and I can send it to you.
library(future.apply)
plan(multisession, workers = 16)
perms <- future_replicate(1e3, fit_permuted_RR(m), FALSE, FALSE)
# Save the permutations, as we will also use them for the other analyses:
write_rds(perms, 'evolrate_RR_perms.rds')

# ratio of evolutionary rates in the permutations:
p_diffs <- map_dbl(perms, ~mean(magnitude(.x$multiple.rates[,4:6]) / mean(magnitude(.x$multiple.rates[,1:3]))))
# observed ratio of evolutionary rates:
o_diff <- mean(magnitude(RR$multiple.rates[,4:6])) / mean(magnitude(RR$multiple.rates[,1:3]))
# p-value:
pmin(mean(p_diffs >= o_diff), mean(p_diffs <= o_diff)) * 2
# observed = 1.26, p = 0.001

## Summarise as figure. S2 in the paper.
ggplot() + 
  geom_vline(aes(xintercept = 1), size = 1, color = 'grey50') +
  geom_density(aes(x, fill = 'simulated'), data.frame(x = p_diffs), alpha = 0.6, col = NA) + 
  geom_vline(aes(xintercept = x, color = 'observed'), data.frame(x = o_diff), size = 1) +
  scale_fill_manual(values = 'grey20') +
  scale_x_log10(
    limits = c(7.5/10, 10/7.5), breaks = c(4/5, 10/9, 1, 9/10, 5/4), 
    labels = expression(frac(4, 5), frac(10, 9), 1, frac(9, 10), frac(5, 4))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = 'firebrick') +
  labs(
    x = expression(frac('mean male rate', 'mean female rate')),
    y = 'density',
    lty = NULL, color = NULL, fill = NULL,
    caption = paste('p =', pmin(mean(o_diff <= p_diffs), mean(o_diff >= p_diffs)) * 2)
  ) +
  theme_minimal() + 
  theme(legend.position = 'top', legend.margin = margin(0, 0, 0, 0))
ggsave('paper_plots/figure_S2.png', dpi = 400, width = 5, height = 2.5)
