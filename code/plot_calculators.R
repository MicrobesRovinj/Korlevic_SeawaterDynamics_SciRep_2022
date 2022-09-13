#################################################################################################################
# plot_calculators.R
# 
# A script to plot richness and diversity calculators.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/calculators.jpg
#           
#################################################################################################################

# Loading OTU/sample data
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared")

# Adding sample labels to row names (input for library vegan)
shared <- shared %>%
  select(-label, -numOtus) %>%
  column_to_rownames("Group")

# Generating a random rarefied community data
rarefied <- shared %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair = "unique", rownames = NA) %>%
  rownames_to_column("Group") %>%
  select_if(list(~ !is.numeric(.) || sum(.) !=0 )) %>%
  column_to_rownames("Group")

# Saving rarefied community data
save(rarefied, file = "results/numerical/rarefied.Rdata")

# Calculating Chao1 and ACE species estimators
estimators <- estimateR(rarefied)
estimators <- as_tibble(t(estimators)) %>%
  add_column("Group" = colnames(estimators), .before = TRUE)

# Calculating diversity indices
shannon <- diversity(rarefied, index = "shannon")
shannon <- tibble(Group = names(shannon), shannon)
invsimpson <- diversity(rarefied, index = "invsimpson")
invsimpson <- tibble(Group = names(invsimpson), invsimpson)

# Transforming Shannon entropy to effective number of species
# (http://www.loujost.com/Statistics%20and%20Physics/Diversity%20and%20Similarity/EffectiveNumberOfSpecies.htm)
eshannon <- mutate(shannon, shannon = exp(shannon)) %>%
  rename(eshannon = shannon)

# Joining together estimators and indices
estimators_indices <- inner_join(estimators, eshannon, by = "Group") %>%
  inner_join(., invsimpson, by = "Group") %>%
  mutate(eshannon = as.numeric(eshannon)) %>%
  mutate(invsimpson = as.numeric(invsimpson))         

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with estimators and indices
Sys.setlocale(locale = "en_GB.utf8")
estimators_indices_metadata <- inner_join(metadata, estimators_indices, by = c("ID" = "Group")) %>%
  mutate(date = as.Date(date, "%d.%m.%Y"))

# Saving calculated estimators and indices
save(estimators_indices_metadata, file = "results/numerical/estimators_indices_metadata.Rdata")

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_line(color = "gray60"),
               axis.ticks = element_line(color = "gray60"),
               axis.text = element_text(size = 12, color = "black"),
               axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 2.5),
               axis.title = element_text(size = 14, color = "black"),
               plot.margin = unit(c(5.5, 16.5, 5.5, 16.5), "pt"),
               legend.position = "none",
               plot.title = element_text(size = 16, hjust = 0.5))

# Tidying data for plotting 
estimators_indices_metadata <- estimators_indices_metadata %>%
  gather("S.obs", "S.chao1", "S.ACE", "eshannon", "invsimpson", key = "estimator_index", value = "value")

# Defining line types, dot shapes and fill dot colours
lines_p1 <- c("S.obs" = "dotted", "S.chao1" = "solid", "S.ACE" = "dotted")
shapes_p1 <- c("S.obs" = 21, "S.chao1" = 23, "S.ACE" = 25)
fills_p1 <- c("S.obs" = "white", "S.chao1" = "black", "S.ACE" = "white")

lines_p2 <- c(c("eshannon" = "solid", "invsimpson" = "dotted"))
shapes_p2 <- c("eshannon" = 21, "invsimpson" = 24)
fills_p2 <- c("eshannon" = "black", "invsimpson" = "white")

# Generating plots
# Saline
saline_p1 <- filter(estimators_indices_metadata, station == "S") %>%
  filter(estimator_index == "S.obs" | estimator_index == "S.chao1" | estimator_index == "S.ACE") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p1) +
  scale_shape_manual(values = shapes_p1) +
  scale_fill_manual(values = fills_p1) +
  scale_y_continuous(limits = c(0, 3500), breaks = c(seq(0, 3500, by = 500)), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "", y = "") +
  ggtitle(parse(text = "bold('Bay of Saline')")) +
  theme +
  theme(axis.text.x = element_blank())

saline_p2 <- filter(estimators_indices_metadata, station == "S") %>%
  filter(estimator_index == "eshannon" | estimator_index == "invsimpson") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p2) +
  scale_shape_manual(values = shapes_p2) +
  scale_fill_manual(values = fills_p2) +
  scale_y_continuous(limits = c(0, 150), breaks = c(seq(0, 150, by = 50)), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "", y = "") +
  theme +
  theme(axis.text.x = element_blank())

# Funtana
funtana_p1 <- filter(estimators_indices_metadata, station == "F") %>%
  filter(estimator_index == "S.obs" | estimator_index == "S.chao1" | estimator_index == "S.ACE") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p1) +
  scale_shape_manual(values = shapes_p1) +
  scale_fill_manual(values = fills_p1) +
  scale_y_continuous(limits = c(0, 3500), breaks = c(seq(0, 3500, by = 500)), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "", y = "") +
  ggtitle(parse(text = "bold('Bay of Funtana')")) +
  theme +
  theme(axis.text.x = element_blank())

funtana_p2 <- filter(estimators_indices_metadata, station == "F") %>%
  filter(estimator_index == "eshannon" | estimator_index == "invsimpson") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p2) +
  scale_shape_manual(values = shapes_p2) +
  scale_fill_manual(values = fills_p2) +
  scale_y_continuous(limits = c(0, 150), breaks = c(seq(0, 150, by = 50)), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "Date", y = "") +
  theme

# Generating a plot to extract a common legend
labels <- c("Observed Number of OTUs", "Chao1", "ACE", "Exponential Shannon", "Inverse Simpson")
p_legend <- filter(estimators_indices_metadata, station == "F") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = c(lines_p1, lines_p2),
                        breaks = c(names(lines_p1), names(lines_p2)),
                        labels = labels) +
  scale_shape_manual(values = c(shapes_p1, shapes_p2),
                     breaks = c(names(lines_p1), names(lines_p2)),
                     labels = labels) +
  scale_fill_manual(values = c(fills_p1, fills_p2),
                    breaks = c(names(lines_p1), names(lines_p2)),
                    labels = labels) +
  labs(x = "", y = "") +
  theme +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(r = 0.2, unit = "cm")),
        legend.key.width = unit(1.4, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.key = element_rect(fill = "white"),
        legend.justification = c("top"),
        legend.text.align = 0) +
  guides(linetype = guide_legend(ncol = 2))
legend <- cowplot :: get_legend(p_legend)

# Combining plots together and saving
p <- cowplot :: plot_grid(saline_p1, saline_p2, funtana_p1, funtana_p2, ncol = 1, nrow = 4,
                          rel_heights = c(0.29, 0.177, 0.29, 0.243), align = "v")
cowplot_theme_1 <- cowplot :: draw_label("Number of OTUs", x = 0.045, y = 0.33,
                                          vjust = -0.5, angle = 90, fontfamily = "Times", size = 14)
cowplot_theme_2 <- cowplot :: draw_label("Number of OTUs", x = 0.045, y = 0.76,
                                          vjust = -0.5, angle = 90, fontfamily = "Times", size = 14)
p <- cowplot :: plot_grid(p, ncol  =  1, nrow  =  2, legend, rel_heights  =  c(4.6, 0.4)) +
  cowplot_theme_1 +
  cowplot_theme_2
ggsave("results/figures/calculators.jpg", p, width = 210, height = 297, units = "mm", bg = "white")
