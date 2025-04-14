### visualize covariance data --> "cov_data.csv" developed manually from evaluation outputs
cov_data <- read.csv("springs_cov_rare_18s.csv") %>%
  filter(cat == "Mount Hood" | cat == "Glacier National Park")

cov_mh <- filter(cov_data, cat == "Mount Hood") # extract MH data
cov_gnp <- filter(cov_data, cat == "Glacier National Park") # extract GNP data

cov_data_order <- cov_data %>% # reorder variables to reflect an increasing mean r_squared
  mutate(env_variable = fct_reorder(env_variable, r_squared, mean, na.rm = TRUE,
                                    .desc = TRUE))
# plot
ggplot(cov_data_order, aes(x = env_variable, y = r_squared, color = cat, shape = cat)) +
  geom_jitter(alpha = 0.7, width = 0.1, size = 5, show.legend = FALSE) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  geom_hline(yintercept = summarise(cov_mh, 
                                    mean(r_squared, na.rm = TRUE)) %>% pull(), 
             lty = 2, color = wes_palette("Moonrise2")[2], lwd = 2)+ 
  geom_hline(yintercept = summarise(cov_gnp, 
                                    mean(r_squared, na.rm = TRUE)) %>% pull(), 
             lty = 2, color = wes_palette("Moonrise2")[1], lwd = 2)+
  labs(y = bquote(~R^2)) +
  stat_summary(fun.data = "mean_cl_boot", #Add mean and 05% confidence intervals
               show.legend = TRUE,
               position = position_nudge(x = .3, y = 0),
               fun.args=(conf.int=0.95),
               size = 1,
               pch = 1) +
  coord_flip() +
  scale_y_continuous(limits = c(0,1)) +
  theme_light() +
  guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  theme(axis.title = element_text(size = 20), #Modify text sizes to improve figure readability
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")

# compare geographic means
cov_data_order %>%
  group_by(cat) %>%
  summarise (mean = mean(r_squared))

# statistical analysis for R2
t.test(r_squared ~ cat, data = cov_data_order)
cov_sum <- cov_data_order %>%
  group_by(env_variable, cat) %>%
  summarise(mean = mean(r_squared))

# range of differences
cov_sum %>%
  group_by(env_variable) %>%
  summarise(diff = diff(mean)) %>%
  filter(env_variable == "Cl" |
           env_variable == "Na" |
           env_variable == "Sr" |
           env_variable == "K" |
           env_variable == "SiO2" |
           env_variable == "SO4" |
           env_variable == "DO" |
           env_variable == "Mg") %>%
  arrange(diff)
