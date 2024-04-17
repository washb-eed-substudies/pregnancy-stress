rm(list=ls())

source(here::here("0-config.R"))

res <- readRDS(here("results/adjusted/adj_res.RDS"))

plotdf <- res %>% mutate(contrast=paste0(Y,"-",X)) %>%
  arrange(coef) %>%
  mutate(contrast=factor(contrast, levels=unique(contrast)))

p_forest <- ggplot(res, aes(x=Y, y=point.diff)) + geom_point() +
  geom_linerange(aes(ymin=lb.diff, ymax=ub.diff)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~X) + coord_flip()

p_forest


