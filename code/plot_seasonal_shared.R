#################################################################################################################
# code/plot_seasonal_shared.R
# 
# A script to plot the Bray-Curtis and Jaccard's Similarity Coefficients of subsequent sampling dates.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/temporal_shared.jpg
#
#################################################################################################################

# Loading OTU/sample data
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared")

# Generating a random rarefied community data
rarefied <- shared %>%
  select(-label, -Group, -numOtus) %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair = "unique") %>%
  add_column("Group" = shared$Group, .before = TRUE) %>%
  select_if(list(~ !is.numeric(.) || sum(.) != 0))

# Copying sample labels to rows (input for library vegan)
rarefied <- rarefied %>%
  column_to_rownames("Group")

# Calculating dissimilarity indices
jaccard <- as.matrix(vegdist(rarefied, method = "jaccard", binary = TRUE))
jaccard <- as_tibble(jaccard) %>%
  add_column("V1" = rownames(jaccard), .before = TRUE) %>%
  gather(key = "V2", value = "jaccard", 2 : ncol(.))

bray <- as.matrix(vegdist(rarefied, method = "bray", binary = FALSE))
bray <- as_tibble(bray) %>%
  add_column("V1" = rownames(bray), .before = TRUE) %>%
  gather(key = "V2", value = "bray", 2 : ncol(.))

distance <- inner_join(jaccard, bray, by = c("V1" = "V1", "V2" = "V2")) %>%
  mutate_at(c("V1", "V2"), list(~ as.character(.)))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with dissimilarity index data
Sys.setlocale(locale = "en_GB.utf8")
distance_metadata <- inner_join(metadata, distance, by = c("ID" = "V1")) %>%
  inner_join(., metadata, by = c("V2" = "ID")) %>%
  mutate(date.x = as.Date(date.x, "%d.%m.%Y")) %>%
  mutate(date.y = as.Date(date.y, "%d.%m.%Y")) %>%
  mutate(jaccard = (1 - jaccard) * 100) %>%
  mutate(bray = (1 - bray) * 100) %>%
  gather(key = "index", value = "value", jaccard, bray)

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

# Defining line types, dot shapes and fill dot colours
lines <- c("jaccard" = "dotted", "bray" = "solid")
shapes <- c("jaccard" = 21, "bray" = 23)
fills <- c("jaccard" = "white", "bray" = "black")

# Generating plots
# Saline
data <- filter(distance_metadata, station.x == "S" & station.y == "S") %>%
  filter((date.x == "2017-07-11" & date.y == "2017-07-26") |
         (date.x == "2017-07-26" & date.y == "2017-08-09") |
         (date.x == "2017-08-09" & date.y == "2017-08-23") |
         (date.x == "2017-08-23" & date.y == "2017-09-14") |
         (date.x == "2017-09-14" & date.y == "2017-10-11") |
         (date.x == "2017-10-11" & date.y == "2017-11-22") |
         (date.x == "2017-11-22" & date.y == "2017-12-13") |
         (date.x == "2017-12-13" & date.y == "2018-02-12") |
         (date.x == "2018-02-12" & date.y == "2018-03-26") |
         (date.x == "2018-03-26" & date.y == "2018-04-23") |
         (date.x == "2018-04-23" & date.y == "2018-05-21") |
         (date.x == "2018-05-21" & date.y == "2018-06-18") |
         (date.x == "2018-06-18" & date.y == "2018-07-09") |
         (date.x == "2018-07-09" & date.y == "2018-08-08") |
         (date.x == "2018-08-08" & date.y == "2018-09-03") |
         (date.x == "2018-09-03" & date.y == "2018-10-04"))

saline_p1 <- filter(data, index == "bray") %>%
  ggplot(aes(x = date.y, y = value, linetype = index, shape = index, fill = index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = fills) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(seq(0, 100, by = 25)), expand = c(0, 0)) +
  labs(x = "", y = "%") +
  ggtitle(parse(text = "bold('Bay of Saline')")) +
  theme +
  theme(axis.text.x = element_blank())

saline_p2 <- filter(data, index == "jaccard") %>%
  ggplot(aes(x = date.y, y = value, linetype = index, shape = index, fill = index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = fills) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                         "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                         "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 30), breaks = c(seq(0, 30, by = 5)), expand = c(0, 0)) +
  labs(x = "", y = "%") +
  theme +
  theme(axis.text.x = element_blank())

# Funtana
data <- filter(distance_metadata, station.x == "F" & station.y == "F") %>%
  filter((date.x == "2017-07-13" & date.y == "2017-07-27") |
         (date.x == "2017-07-27" & date.y == "2017-08-10") |
         (date.x == "2017-08-10" & date.y == "2017-08-24") |
         (date.x == "2017-08-24" & date.y == "2017-09-19") |
         (date.x == "2017-09-19" & date.y == "2017-10-12") |
         (date.x == "2017-10-12" & date.y == "2017-11-23") |
         (date.x == "2017-11-23" & date.y == "2017-12-14") |
         (date.x == "2017-12-14" & date.y == "2018-02-13") |
         (date.x == "2018-02-13" & date.y == "2018-03-27") |
         (date.x == "2018-03-27" & date.y == "2018-04-24") |
         (date.x == "2018-04-24" & date.y == "2018-05-22") |
         (date.x == "2018-05-22" & date.y == "2018-06-19") |
         (date.x == "2018-06-19" & date.y == "2018-07-10") |
         (date.x == "2018-07-10" & date.y == "2018-08-09") |
         (date.x == "2018-08-09" & date.y == "2018-09-04") |
         (date.x == "2018-09-04" & date.y == "2018-10-05"))

funtana_p1 <- filter(data, index == "bray") %>%
  ggplot(aes(x = date.y, y = value, linetype = index, shape = index, fill = index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = fills) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(seq(0, 100, by = 25)), expand = c(0, 0)) +
  labs(x = "", y = "%") +
  ggtitle(parse(text = "bold('Bay of Funtana')")) +
  theme +
  theme(axis.text.x = element_blank())

funtana_p2 <- filter(data, index == "jaccard") %>%
  ggplot(aes(x = date.y, y = value, linetype = index, shape = index, fill = index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = fills) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 30), breaks = c(seq(0, 30, by = 5)), expand = c(0, 0)) +
  labs(x = "Date", y = "%") +
  theme

# Generating a plot to extract a common legend
labels <- c("Bray-Curtis Similarity Coefficient", "Jaccard's Similarity Coefficient")
p_legend <- filter(distance_metadata, station.x == "F" & station.y == "F") %>%
  ggplot(aes(x = date.y, y = value, linetype = index, shape = index, fill = index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines, labels = labels) +
  scale_shape_manual(values = shapes, labels = labels) +
  scale_fill_manual(values = fills, labels = labels) +
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
  guides(linetype = guide_legend(ncol = 1))
legend <- cowplot :: get_legend(p_legend)

# Combining plots together and saving
p <- cowplot :: plot_grid(saline_p1, saline_p2, funtana_p1, funtana_p2, ncol = 1, nrow = 4,
                        rel_heights = c(0.29, 0.177, 0.29, 0.243), align = "v")
p <- cowplot :: plot_grid(p, ncol = 1, nrow = 2, legend, rel_heights = c(4.6, 0.4))
ggsave("results/figures/seasonal_shared.jpg", p, width = 210, height = 297, units = "mm", bg = "white")
