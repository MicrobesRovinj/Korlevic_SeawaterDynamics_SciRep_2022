#################################################################################################################
# plot_rarefaction.R
# 
# A script to plot the rarefaction curve of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction
#               data/raw/metadata.csv
# Produces: results/figures/rarefaction.jpg
#
#################################################################################################################

# Loading input data and selection of values for plotting
rarefaction <- read_tsv(file = "data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  gather(-numsampled, key = sample, value = sobs) %>%
  mutate(sample = str_replace_all(sample, pattern = "0.03-", replacement = "")) %>%
  drop_na()

# Loading metadata
Sys.setlocale(locale = "en_GB.utf8")
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  arrange(date) %>%
  mutate(date = format(date, "%d %B %Y")) %>%
  mutate(date = str_replace(date, "^0", "")) %>%
  mutate(date = if_else(label == "24/4-18 F-1", paste0(date, "a"), date)) %>%
  mutate(date = if_else(label == "24/4-18 F-2", paste0(date, "b"), date))

# Joining metadata and input data
metadata_rarefaction <- inner_join(metadata, rarefaction, by = c("ID" = "sample"))

# Defining line colour and type for each sampling date
color_type <- tribble(~date, ~color, ~line_type,
                      "11 July 2017", "#A6CEE3", 1,
                      "13 July 2017", "#A6CEE3", 1,
                      "26 July 2017", "#1F78B4", 1,
                      "27 July 2017", "#1F78B4", 1,
                      "9 August 2017", "#B2DF8A", 1,
                      "10 August 2017", "#B2DF8A", 1,
                      "23 August 2017", "#33A02C", 1,
                      "24 August 2017", "#33A02C", 1,
                      "14 September 2017", "#FB9A99", 1,
                      "19 September 2017", "#FB9A99", 1,
                      "11 October 2017", "#E31A1C", 1,
                      "12 October 2017", "#E31A1C", 1,
                      "22 November 2017", "#FDBF6F", 1,
                      "23 November 2017", "#FDBF6F", 1,
                      "13 December 2017", "#CAB2D6", 1,
                      "14 December 2017", "#CAB2D6", 1,
                      "12 February 2018", "#6A3D9A", 1,
                      "13 February 2018", "#6A3D9A", 1,
                      "26 March 2018", "#A6CEE3", 3,
                      "27 March 2018", "#A6CEE3", 3,
                      "23 April 2018", "#1F78B4", 3,
                      "24 April 2018a", "#1F78B4", 3,
                      "24 April 2018b", "#B2DF8A", 3,
                      "21 May 2018", "#33A02C", 3,
                      "22 May 2018", "#33A02C", 3,
                      "18 June 2018", "#FB9A99", 3,
                      "19 June 2018", "#FB9A99", 3,
                      "9 July 2018", "#E31A1C", 3,
                      "10 July 2018", "#E31A1C", 3,
                      "8 August 2018", "#FDBF6F", 3,
                      "9 August 2018", "#FDBF6F", 3,
                      "3 September 2018", "#FF7F00", 3,
                      "4 September 2018", "#FF7F00", 3,
                      "4 October 2018", "#CAB2D6", 3,
                      "5 October 2018", "#CAB2D6", 3)

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_line(color = "gray60"),
               axis.ticks = element_line(color = "gray60"),
               axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 14, color = "black"),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               legend.position = "right",
               legend.title = element_blank(),
               legend.text = element_text(size = 10, margin = margin(r = 0.2, unit = "cm")),
               legend.key.width = unit(1.4, "cm"),
               legend.key.height = unit(0.5, "cm"),
               legend.key = element_rect(fill = "white"),
               legend.justification = c("top"),
               legend.text.align = 0,
               legend.spacing.x = unit(0, "cm"),
               plot.title = element_text(size = 16, hjust = 0.5))

# Plots generation
p1 <- filter(metadata_rarefaction, station == "S") %>%
  ggplot(aes(x = numsampled, y = sobs, group = ID, color = factor(date, levels = unique(date)),
             linetype = factor(date, levels = unique(date)))) +
  geom_line(size = 1.5) +
  scale_colour_manual(name = NULL,
                      breaks = color_type$date,
                      values = color_type$color) +
  scale_linetype_manual(name = NULL,
                        breaks = color_type$date,
                        values = color_type$line_type) +
  scale_x_continuous(limits = c(0, 70000), breaks = c(seq(0, 70000, by = 10000)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2000), breaks = c(seq(0, 2000, by = 500)), expand = c(0, 0)) +
  labs(x = "", y = "Number of OTUs") +
  ggtitle(parse(text = "bold('Bay of Saline')")) +
  theme

p2 <- filter(metadata_rarefaction, station == "F") %>%
  ggplot(aes(x = numsampled, y = sobs, group = ID, color = factor(date, levels = unique(date)),
             linetype = factor(date, levels = unique(date)))) +
  geom_line(size = 1.5) +
  scale_colour_manual(name = NULL,
                      breaks = color_type$date,
                      values = color_type$color) +
  scale_linetype_manual(name = NULL,
                        breaks = color_type$date,
                        values = color_type$line_type) +
  scale_x_continuous(limits = c(0, 80000), breaks = c(seq(0, 80000, by = 10000)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2500), breaks = c(seq(0, 2500, by = 500)), expand = c(0, 0)) +
  labs(x = "Number of Sequences", y = "Number of OTUs") +
  ggtitle(parse(text = "bold('Bay of Funtana')")) +
  theme

# Combining plots together and save
p <- cowplot :: plot_grid(p1, p2, ncol = 1, nrow = 2)
ggsave("results/figures/rarefaction.jpg", p, width = 210, height = 297, units = "mm")
