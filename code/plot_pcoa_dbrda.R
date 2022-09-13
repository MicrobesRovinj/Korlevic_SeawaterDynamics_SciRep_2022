#################################################################################################################
# plot_pcoa_dbrda.R
# 
# A script to generate the PCoA and db-RDA figures.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/environmental_data.csv
#               data/raw/metadata.csv
# Produces: results/figures/pcoa_dbrda_figure.jpg
#
#################################################################################################################

#######################
# PCoA
#######################

# Loading rarefied community data created during the generation of the richness and diversity calculators plot
load(file = "results/numerical/rarefied.Rdata")
rarefied <- rarefied %>%
  rownames_to_column("Group")

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Formatting sampling months
metadata <- metadata %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(month = format(date, "%B")) %>%
  mutate(year = format(date, "%Y"))

# Formatting sampling stations
metadata <- metadata %>%
  mutate(station = str_replace(station, "^S$", "Saline")) %>%
  mutate(station = str_replace(station, "^F$", "Funtana"))

# Joining metadata and OTU/sample data
rarefied_metadata <- inner_join(metadata, rarefied, by = c("ID" = "Group"))

# Function for setting the number of decimal places
scaleFUN <- function(x) sprintf("%.1f", x)

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 18, color = "black"),
               plot.margin = unit(c(5.5, 16.5, 16.5, 16.5), "pt"),
               legend.text = element_text(size = 18, margin = margin(r = 0.2, unit = "cm")),
               legend.text.align = 0,
               legend.key = element_rect(fill = NA),
               legend.key.width = unit(0, "cm"),
               legend.key.height = unit(1.0, "cm"),
               legend.margin = margin(-5, 0, 0, 0),
               legend.justification = c("left", "bottom"))

# Defining colours for each month
colours_months <- c("February" = "#A6CEE3",
                    "March" = "#1F78B4",
                    "April" = "#B2DF8A",
                    "May" = "#33A02C",
                    "June" = "#FB9A99",
                    "July" = "#E31A1C" ,
                    "August" = "#FDBF6F",
                    "September" = "#FF7F00",
                    "October" =  "#CAB2D6",
                    "November" = "#6A3D9A",
                    "December" = "#B15928")

# Defining shapes for each station and year
shapes_stations_years <- c("2017" = 3,
                           "2018" = NA,
                           "Saline" = 21,
                           "Funtana" = 23)

# Selecting samples for plotting
rarefied_metadata_select <- rarefied_metadata

# Generating PCoA data
spe.bray <- rarefied_metadata_select %>%
  column_to_rownames("ID") %>%
  select(starts_with("Otu")) %>%
  select(where(~ sum(.) !=  0)) %>%
  vegdist(method = "bray")

spe.b.pcoa <- wcmdscale(spe.bray, k = (nrow(rarefied_metadata_select) - 1), eig = TRUE, add = "lingoes")

coordinates <- spe.b.pcoa %>%
  scores(choices = c(1, 2)) %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  rename(A1 = Dim1, A2 = Dim2) %>%
  as_tibble()

coordinates <- inner_join(metadata, coordinates, by = c("ID" = "ID"))

# Generating the PCoA plot
p.pcoa <- ggplot() +
  geom_point(data = coordinates, aes(x = A1, y = A2, fill = month, shape = station), size = 5, stroke = 0.5) +
  geom_point(data = coordinates, aes(x = A1, y = A2, shape = year), size = 7, stroke = 0.3) +
  scale_fill_manual(name = NULL,
                    breaks = names(colours_months),
                    values = colours_months) +
  scale_shape_manual(name = NULL,
                    breaks = names(shapes_stations_years),
                    values = shapes_stations_years) +
  labs(x = paste0("PCoA I (", format(round(spe.b.pcoa$eig[1] / sum(spe.b.pcoa$eig) * 100, digits = 1), nsmall = 1), " %)"), 
       y = paste0("PCoA II (", format(round(spe.b.pcoa$eig[2] / sum(spe.b.pcoa$eig) * 100, digits = 1), nsmall = 1), " %)")) +
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  coord_fixed(ratio = 1, xlim = c(-0.30, 0.50), ylim = c(-0.30, 0.50)) +
  theme +
  theme(legend.position = "")

#######################
# db-RDA
#######################

# The DNA sample originating from Funtana sampled on 24 April 2018 (ID, 23) was sequences twice.
# Sequences from both samples are summed.
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared") %>%
  mutate(Group = str_replace(Group, "23_.", "23")) %>%
  group_by(Group) %>%
  summarise(across(starts_with("Otu"), list(sum), .names = "{.col}"), .groups = "drop")

# Generating the random rarefied community data
rarefied <- shared %>%
  column_to_rownames("Group") %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair = "unique", rownames = NA) %>%
  rownames_to_column("Group") %>%
  select(!where(is.numeric) | where(~ is.numeric(.) && sum(.) != 0))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(ID = str_replace(ID, "23_1", "23"))

# Formatting sampling months
metadata <- metadata %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(month = format(date, "%B")) %>%
  mutate(year = format(date, "%Y"))

# Formatting sampling stations
metadata <- metadata %>%
  mutate(station = str_replace(station, "^S$", "Saline")) %>%
  mutate(station = str_replace(station, "^F$", "Funtana"))

# Joining metadata and OTU/sample data
rarefied_metadata <- inner_join(metadata, rarefied, by = c("ID" = "Group"))

# Loading environmental variables data
env <- read_tsv("data/raw/environmental_data.csv") %>%
  mutate(ID = as.character(ID)) %>%
  rename_all(~ str_replace_all(., " \\(.*\\)", ""))

# Joining metadata and environmental variables data
env_metadata <- inner_join(metadata, env, by = c("ID" = "ID"))

# Selecting samples for plotting
rarefied_metadata_select <- rarefied_metadata

# Generating db-RDA data
rarefied_metadata_select <- rarefied_metadata %>%
  column_to_rownames("ID") %>%
  select(starts_with("Otu")) %>%
  select(where(~ sum(.) !=  0))

spe.dbrda <- capscale(rarefied_metadata_select ~ `T` + S + PO4 + NH4 + NO2 + NO3 + SiO4 + PM + Chla + PA, data = env_metadata,
                      distance = "bray", add = "lingoes", dfun = vegdist, na.action = na.fail)

(R2adj <- RsquareAdj(spe.dbrda)$adj.r.squared)

# Global test of the db-RDA result
(anova.global <- anova.cca(spe.dbrda, permuatations = how(nperm = 999)))

# Test of all canonical axes
(anova.axes <- anova(spe.dbrda, by = "axis", permutations = how(nperm = 999)))

# Variance inflation factors (VIF)
vif.cca(spe.dbrda)

# Combining and saving statistic results
statistic.dbrda <- list(R2adj = R2adj, anova.global = anova.global, anova.axes = anova.axes)
save(statistic.dbrda, file = "results/numerical/statistic.dbrda.Rdata")

# Extracting data to generate the db-RDA plot
coordinates_sites <- spe.dbrda %>%
  scores(choices = c(1, 2), display = "lc", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  rename(A1 = CAP1, A2 = CAP2) %>%
  as_tibble()
coordinates_sites <- inner_join(metadata, coordinates_sites, by = c("ID" = "ID"))

coordinates_variables <- spe.dbrda %>%
  scores(choices = c(1, 2), display = "bp", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("variables") %>%
  rename(A1 = CAP1, A2 = CAP2) %>%
  as_tibble() %>%
  mutate(variables = case_when(variables == "PO4" ~ "bold('PO'['4']^{' 3–'})",
                               variables == "NH4" ~ "bold('NH'['4']^{' +'})",
                               variables == "NO2" ~ "bold('NO'['2']^{' –'})",
                               variables == "NO3" ~ "bold('NO'['3']^{' –'})",
                               variables == "SiO4" ~ "bold('Si(OH)'['4'])",
                               variables == "Chla" ~ "bold('Chl')~bolditalic('a')",
                               TRUE ~ paste0("bold(", variables, ")")))

eigenvalues_constrained <- spe.dbrda %>%
  summary(scaling = 2) %>%
  .$concont %>%
  .$importance %>%
  as.data.frame() %>%
  rename(A1 = CAP1, A2 = CAP2)

# Generating the db-RDA plot
p.dbrda <- ggplot() +
  geom_point(data = coordinates_sites, aes(x = A1, y = A2, fill = month, shape = station), size = 5, stroke = 0.5) +
  geom_point(data = coordinates_sites, aes(x = A1, y = A2, shape = year), size = 7, stroke = 0.3) +
  scale_fill_manual(name = NULL,
                    breaks = names(colours_months),
                    values = colours_months) +
  scale_shape_manual(name = NULL,
                     breaks = names(shapes_stations_years),
                     values = shapes_stations_years) +
  geom_segment(data = coordinates_variables, aes(x = 0, y = 0, xend = A1, yend = A2),
               arrow = arrow(length = unit(0.25, "cm")),  size = 0.5, color = "gray50") +
  geom_text(data = coordinates_variables,
            aes(x = A1 + 0.15 * cos(atan(A2 / A1)) * sign(A1),
                y = A2 + 0.15 * sin(atan(A2 / A1)) * sign(A1),
                label = variables, family = "Times"),
            col = "gray50", size = 5, parse = TRUE) +
  annotate("text", x = 1.11, y = 1.45,
           label = paste0("bolditalic('R'['a']^{bold('2')})~bold('=')~",
                 "bold('", format(round(R2adj*100, digits = 1), nsmall = 1), "')", "~bold('%')"),
           family = "Times", parse = TRUE, col = "gray50", size = 9) + 
  labs(x = paste0("db-RDA I (", format(round(eigenvalues_constrained$A1[2] * R2adj * 100, digits = 1), nsmall = 1), " %)"), 
       y = paste0("db-RDA II (", format(round(eigenvalues_constrained$A2[2] * R2adj * 100, digits = 1), nsmall = 1), " %)")) +
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  coord_fixed(ratio = 1, xlim = c(-1.10, 1.50), ylim = c(-1.10, 1.50)) +
  theme +
  theme(legend.position = "")

# Extracting a common legend
p <- p.dbrda +
  theme(legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  guides(fill = guide_legend(override.aes = list(shape = 21), order = 1), shape = guide_legend(override.aes = list(size = 5), order = 2))
legend <- cowplot :: get_legend(p)

# Combining plots and saving
p <- cowplot :: plot_grid(p.pcoa, p.dbrda, labels = c("a", "b"), label_size = 48, label_fontfamily = "Times",
                        hjust = -0.1, vjust = 2.0)
p <- cowplot :: ggdraw() +
  cowplot :: draw_plot(p, x = 0, y = 0, width = 0.92, height = 1) +
  cowplot :: draw_plot(legend, x = 0.92, y = 0.12, width = 0, height = 1)
ggsave("results/figures/pcoa_dbrda_figure.jpg", p, width = 1.5 * 297, height = 210, units = "mm", bg = "white")
