## This script was used to compare the color structure of the field guide drawings to photographs
## of museum specimens

library(tidyverse)
short_list <- data.table::fread('centroids_photos.csv')

# All pairwise comparison
dE_photo_all <- short_list %>% 
  mutate(Tolman = paste(Tolman, sex, sep = '_')) %>% 
  select(Tolman, L, a, b) %>% 
  expand_grid(a = ., b = .) %>% 
  flatten() %>% as.data.frame() %>% 
  setNames(c('Tolman', 'L', 'a', 'b', 'other_Tolman', 'L2', 'a2', 'b2')) %>% 
  as_tibble() %>% 
  mutate(
    L_photo = L - L2, a_photo = a - a2, b_photo = b - b2
  ) %>% 
  select(Tolman, other_Tolman, ends_with('_photo'))

dE_drawing_all <- data.table::fread('centroids.csv') %>% 
  select(Tolman = Species, L_mean_females, L_mean_males, a_mean_females, a_mean_males, b_mean_females, b_mean_males) %>% 
  gather(var, val, -Tolman) %>% 
  mutate(var = str_remove(var, '_mean')) %>% 
  separate(var, into = c('var', 'sex')) %>% 
  spread(var, val) %>% 
  filter(Tolman %in% short_list$Tolman) %>% 
  mutate(
    sex = str_remove(sex, 's'),
    Tolman = paste(Tolman, sex, sep = '_')
  ) %>% 
  select(Tolman, L, a, b) %>% 
  expand_grid(a = ., b = .) %>% 
  flatten() %>% as.data.frame() %>% 
  setNames(c('Tolman', 'L', 'a', 'b', 'other_Tolman', 'L2', 'a2', 'b2')) %>% 
  as_tibble() %>% 
  mutate(
    L_drawing = L - L2, a_drawing = a - a2, b_drawing = b - b2
  ) %>% 
  select(Tolman, other_Tolman, ends_with('drawing'))

dE_all <- full_join(dE_photo_all, dE_drawing_all, by = c('Tolman', 'other_Tolman')) %>% 
  filter(as.numeric(factor(Tolman)) > as.numeric(factor(other_Tolman)))
dE_all2 <- dE_all %>% 
  gather(key, val, -Tolman, -other_Tolman) %>% 
  separate(key, c('axis', 'method')) %>% 
  spread(method, val) %>% 
  mutate(axis = factor(axis, c('L', 'a', 'b')))

dE_all2 %>% 
  group_by(axis) %>% 
  summarise(r = cor(drawing, photo, method = 'pearson'))
# # A tibble: 3 x 2
#   axis  r
#   <fct> <dbl>
# 1 L     0.873   # also has largest variance
# 2 a     0.849
# 3 b     0.778
dE_all2 %>% 
  filter(axis != 'dE') %>% 
  group_by(axis) %>% 
  summarise(slope = coef(lm(photo ~ drawing))[2])
# # A tibble: 3 x 2
#   axis  slope
#   <fct> <dbl>
# 1 L     0.549
# 2 a     0.435
# 3 b     0.568

# Make Figure S1
split(dE_all2, dE_all2$axis)[-4] %>% 
  map(~ggplot(., aes(drawing, photo)) +
        geom_point(shape = '.') + #geom_smooth(method = 'lm') +
        labs(
          x = 'measured pairwise difference between specimens\n(drawing)', 
          y = 'measured pairwise difference between specimens\n(photo)'
        ) +
        theme_minimal() +
        xlim(min(c(.$drawing, .$photo)), max(c(.$drawing, .$photo))) +
        ylim(min(c(.$drawing, .$photo)), max(c(.$drawing, .$photo)))
  ) %>% 
  map2(map(c('L, r = 0.87', 'a, r = 0.85', 'b, r = 0.78'), ggtitle), `+`) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 3)
ggsave('paper_plots/color_verification.png', width = 4, height = 12)
