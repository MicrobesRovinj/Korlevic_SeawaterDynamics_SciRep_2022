#################################################################################################################
# plot_community_bar_plot_month.R
#
# A script to plot the community structure of each group of months.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/community_bar_plot_month.jpg
#           results/figures/community_bar_plot_month_taxa.jpg
#
#################################################################################################################

# Loading input data containing sequence abundances and subsequent input data customisation
community <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary") %>%
  filter(!str_detect(taxon, "^Eukaryota")) %>%
  filter(taxon != "Root")

# Removing chloroplast and mitochondrial sequences and subtracting their number from higher taxonomic levels to which they belong
chloroplast <- filter(community, str_detect(taxon, "^Chloroplast$"))$rankID
mitochondria <- filter(community, str_detect(taxon, "^Mitochondria$"))$rankID
community <- mutate_at(community, 5 : ncol(community), list(~ case_when(
  rankID == str_extract(chloroplast, "(\\d+\\.){3}\\d+") ~ . - .[taxon == "Chloroplast"],
  rankID == str_extract(chloroplast, "(\\d+\\.){2}\\d+") ~ . - .[taxon == "Chloroplast"],
  rankID == str_extract(chloroplast, "(\\d+\\.){1}\\d+") ~ . - .[taxon == "Chloroplast"],
  TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Chloroplast")) %>%
  mutate_at(5 : ncol(.), list(~ case_when(
    rankID == str_extract(mitochondria, "(\\d+\\.){4}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){3}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){2}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){1}\\d+") ~ . - .[taxon == "Mitochondria"],
    TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Mitochondria")) %>%
  select(-ATCC_1, -ATCC_5, -NC_1) %>%
  mutate(`23` = `23_1` + `23_2`) %>%
  select(-`23_1`, -`23_2`) %>%
  group_by(taxlevel) %>%
  mutate_at(5 : ncol(.), list(~ . / sum(.) * 100)) %>%
  ungroup()

# Replacing the name Campilobacterota with Campylobacterota
community <- community %>%
  mutate(taxon = if_else(taxon == "Campilobacterota", "Campylobacterota", taxon))

# Loading colours for each taxa on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(- Taxlevel) %>%
  deframe()

# Metadata loading and customisation 
metadata <- read_tsv("data/raw/metadata.csv") %>%
  filter(ID != "23_1") %>%
  mutate(ID = str_replace(ID, "23_2", "23")) %>%
  mutate(label = str_replace(label, "24/4-18 F-2", "24/4-18 F")) %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(period = case_when((date >= "2017-06-01" & date < "2017-10-31") |
                            (date >= "2018-06-01" & date < "2018-10-31")~ "June – October",
                            date >= "2017-11-01" & date < "2018-03-31" ~ "November – March",
                            date >= "2018-04-01" & date < "2018-05-31" ~ "April – May",
                            TRUE ~ "")) %>%
  mutate(period = factor(period, levels = c("June – October", "November – March", "April – May")))

# Defining plot titles
title_labs <- c("total_community" = "bold('Total Community')",
                "Alphaproteobacteria" = "bolditalic('Alphaproteobacteria')",
                "Bacteroidota" = "bolditalic('Bacteroidota')",
                "Gammaproteobacteria" = "bolditalic('Gammaproteobacteria')",
                "Cyanobacteria" = "bolditalic('Cyanobacteria')")

#################################################################################################################
# Generating input data for taxonomy bar plots in a for loop
#################################################################################################################

# Defining groups, their corresponding taxonomic levels and threshold values
groups <- tribble(~ taxa, ~ taxlevel, ~ threshold_value,
                  "total_community", 1, 1,
                  "Alphaproteobacteria", 3, 2,
                  "Bacteroidota", 2, 1,
                  "Gammaproteobacteria", 3, 1,
                  "Cyanobacteria", 2, 1)

for (i in groups$taxa) {
  
  if (i == "total_community") {
    
    # Selecting taxonomic groups for plotting
    plot <- filter(community,
                   taxlevel == 2 |
                   (taxlevel == 3 & str_detect(rankID, paste0(filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, "\\.")))) %>%
      filter_at(6 : ncol(.), any_vars(. >= filter(groups, taxa == i)$threshold_value)) %>%
      mutate_at(5 : ncol(.), list(~ case_when(taxon == "Proteobacteria" ~ . - sum(.[taxlevel == 3 & str_detect(rankID, paste0(filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, "\\."))]), TRUE ~ .))) %>%
      mutate(taxon = str_replace(taxon, "Proteobacteria", "Other_Proteobacteria")) %>%
      mutate(taxon = str_replace_all(taxon, c("unknown_unclassified" = "No_Relative", "unknown" = "No_Relative"))) %>%
      filter_at(6 : ncol(.), any_vars(. >= filter(groups, taxa == i)$threshold_value)) %>%
      bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), "Other")))) %>%
      arrange(taxon %in% "No_Relative")
    
    } else {
      
      # Selecting taxonomic groups above threshold value
      select <- filter(community,
                       taxlevel == 6 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, paste0("^", i, "$")))$rankID, "\\."))) %>%
        filter(if_any(6 : ncol(.), ~ . >= filter(groups, taxa == i)$threshold_value))
      
      # Selecting taxonomic groups for plotting
      plot <- filter(community,
                     taxlevel == 6 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, paste0("^", i, "$")))$rankID, "\\."))) %>%
        mutate(across(5 : ncol(.), ~ . / sum(.) * 100)) %>%
        filter(rankID %in% select$rankID) %>%
        bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), paste0("Other_", i))))) %>%
        # Removing last digit from "uncultured" rankID to be used in the next step
        mutate(rankID = if_else(taxon == "uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID)) %>%
        # Removing last two digits from "uncultured_fa" rankID to be used in the next step
        mutate(rankID = if_else(taxon == "uncultured_fa", str_replace(rankID, "(\\.\\d+){2}$", ""), rankID)) %>%
        # Removing last digit from "Unknown_Family" rankID to be used in the next step
        mutate(rankID = if_else(taxon == "Unknown_Family", str_replace(rankID, "\\.\\d+$", ""), rankID))
      
      # Adding information to "uncultured" taxa describing higher taxonomic level to which they belong
      uncultured <- select(filter(community, rankID %in% filter(plot, taxon == "uncultured" |
                                                                      taxon == "uncultured_fa" |
                                                                      taxon == "Unknown_Family")$rankID), rankID, taxon) %>%
        rename(rankID_uncultured = rankID, taxon_uncultured = taxon)
      plot <- left_join(plot, uncultured, by = c("rankID" = "rankID_uncultured")) %>%
        mutate(taxon = if_else(taxon == "uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
        mutate(taxon = if_else(taxon == "uncultured_fa", paste0("uncultured", "_", taxon_uncultured), taxon)) %>%
        mutate(taxon = if_else(taxon == "Unknown_Family", taxon_uncultured, taxon)) %>%
        select(-taxon_uncultured)
      
    }
  
  # Customisation of taxonomic names in legend
  names <- case_when(i == "total_community" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                     str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                plot$taxon == "Bacteria_unclassified" ~ "italic('Bacteria')~plain('(NR)')",
                                                                                plot$taxon == "Marinimicrobia_(SAR406_clade)" ~ "italic('Marinimicrobia')",
                                                                                plot$taxon == "Other" ~ paste0("plain('", plot$taxon, "')"),
                                                                                plot$taxon == "No_Relative" ~ "plain('No Relative')",
                                                                                TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Alphaproteobacteria" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                         str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                    plot$taxon == "OCS116_clade_ge" ~ "plain('OCS116 Clade')",
                                                                                    plot$taxon == "Candidatus_Puniceispirillum" ~ "plain('\"')*italic('Candidatus')~plain('Puniceispirillum\"')",
                                                                                    plot$taxon == "SAR116_clade_ge" ~ "plain('SAR116 Clade')",
                                                                                    plot$taxon == "Stappiaceae_ge" ~ "italic('Stappiaceae')",
                                                                                    plot$taxon == "HIMB11" ~ "plain('HIMB11')",
                                                                                    plot$taxon == "AEGEAN-169_marine_group_ge" ~ "plain('AEGEAN-169 Marine Group')",
                                                                                    plot$taxon == "Clade_Ia" ~ "plain('Subclade Ia')",
                                                                                    plot$taxon == "Clade_Ib" ~ "plain('Subclade Ib')",
                                                                                    plot$taxon == "Clade_II_ge" ~ "plain('Subclade II')",
                                                                                    plot$taxon == "Clade_III_ge" ~ "plain('Subclade III')",
                                                                                    plot$taxon == "Other_Alphaproteobacteria" ~ "plain('Other')~italic('Alphaproteobacteria')",
                                                                                    plot$taxon == "S25-593_ge" ~ "plain('S25-593')",
                                                                                    TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Bacteroidota" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                  str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                             plot$taxon == "NS4_marine_group" ~ "plain('NS4 Marine Group')",
                                                                             plot$taxon == "NS5_marine_group" ~ "plain('NS5 Marine Group')",
                                                                             plot$taxon == "NS11-12_marine_group_ge" ~ "plain('NS11-12 Marine Group')",
                                                                             plot$taxon == "Flavobacteriaceae_unclassified" ~ "italic('Flavobacteriaceae')~plain('(NR)')",
                                                                             plot$taxon == "NS9_marine_group_ge" ~ "plain('NS9 Marine Group')",
                                                                             plot$taxon == "Other_Bacteroidota" ~ "plain('Other')~italic('Bacteroidota')",
                                                                             TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Gammaproteobacteria" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                 str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                            plot$taxon == "OM43_clade" ~ "plain('OM43 Clade')",
                                                                            plot$taxon == "OM60(NOR5)_clade" ~ "plain('OM60 (NOR5) Clade')",
                                                                            plot$taxon == "SAR86_clade_ge" ~ "plain('SAR86 Clade')",
                                                                            plot$taxon == "SUP05_cluster" ~ "plain('SUP05 Cluster')",
                                                                            plot$taxon == "SAR92_clade" ~ "plain('SAR92 Clade')",
                                                                            plot$taxon == "KI89A_clade_ge" ~ "plain('KI89A Clade')",
                                                                            plot$taxon == "OM182_clade_ge" ~ "plain('OM182 Clade')",
                                                                            plot$taxon == "Other_Gammaproteobacteria" ~ "plain('Other')~italic('Gammaproteobacteria')",
                                                                            TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Cyanobacteria" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                     str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                plot$taxon == "Spirulina_DRTO-55.2" ~ "italic('Spirulina')",
                                                                                plot$taxon == "Pleurocapsa_PCC-7319" ~ "italic('Pleurocapsa')",
                                                                                plot$taxon == "Cyanobacteriia_unclassified" ~ "italic('Cyanobacteriia')~plain('(NR)')",
                                                                                plot$taxon == "Nodosilineaceae_unclassified" ~ "italic('Nodosilineaceae')~plain('(NR)')",
                                                                                plot$taxon == "Acrophormium_PCC-7375" ~ "italic('Acrophormium')",
                                                                                plot$taxon == "Phormidesmis_ANT.LACV5.1" ~ "italic('Phormidesmis')",
                                                                                plot$taxon == "Phormidium_MBIC10003" ~ "italic('Phormidium')",
                                                                                plot$taxon == "Cyanobium_PCC-6307" ~ "italic('Cyanobium')",
                                                                                plot$taxon == "Synechococcus_CC9902" ~ "italic('Synechococcus')",
                                                                                plot$taxon == "Schizothrix_LEGE_07164" ~ "italic('Schizothrix')",
                                                                                plot$taxon == "Thermosynechococcales_unclassified" ~ "italic('Thermosynechococcales')~plain('(NR)')",  
                                                                                plot$taxon == "Other_Cyanobacteria" ~ "plain('Other')~italic('Cyanobacteria')",
                                                                                TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     TRUE ~ parse(text = paste0("italic('", plot$taxon, "')")))

  # Tidying sequence abundance data
  plot <- plot %>%
    gather(key = "Group", value = "abundance", 6 : ncol(.))
  
  # Joining sequence abundance data and metadata
  plot <- inner_join(metadata, plot, by = c("ID" = "Group")) %>%
    mutate(taxon = factor(taxon, levels = unique(plot$taxon)))

  if (i == "total_community") {
    
    # Calculating relative abundance of all groups (100 %) in the whole community, tidying the obtained data and joining with metadata
    whole <- community %>%
      filter(taxlevel == filter(groups, taxa == i)$taxlevel) %>%
      bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), sum(.), "Total")))) %>%
      gather(key = "Group", value = "abundance", 6 : ncol(.)) %>%
      filter(taxon == "Total")
    whole <- inner_join(metadata, whole, by = c("ID" = "Group"))
  
  } else {
    
    # Selecting relative abundance of target group in the whole community, tidying the obtained data and joining with metadata
    whole <- community %>%
      filter(taxlevel == filter(groups, taxa == i)$taxlevel) %>%
      gather(key = "Group", value = "abundance", 6 : ncol(.)) %>%
      filter(taxon == i)
    whole <- inner_join(metadata, whole, by = c("ID" = "Group"))
    
  }

  # Generating mean abundance for different groups of months
  plot <- plot %>%
    group_by(period, taxon) %>%
    summarise(mean = mean(abundance), sd = sd(abundance), .groups = "drop")
  
  # Generating mean abundance in the whole community for different groups of months
  whole <- whole %>%
    group_by(period, taxon) %>%
    summarise(mean = mean(abundance), sd = sd(abundance), .groups = "drop")

  # Creating a variable for faceting
  plot <- plot %>%
    mutate(facet = i)
  whole <- whole %>%
    mutate(facet = i)

  # Setting names for objects created in loop
  assign(paste0("plot", sep = "_", i), plot)
  assign(paste0("whole", sep = "_", i), whole)
  assign(paste0("names", sep = "_", i), names)
  
}

# Combining and saving tidied taxonomic results
tidied_taxonomy <- full_join(plot_total_community, plot_Alphaproteobacteria) %>%
  full_join(., plot_Bacteroidota) %>%
  full_join(., plot_Gammaproteobacteria) %>%
  full_join(., plot_Cyanobacteria)
save(tidied_taxonomy, file = "results/numerical/tidied_taxonomy.Rdata")

#################################################################################################################
# Generating taxonomy bar plots
#################################################################################################################

# Generating theme for ggplot
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.background = element_blank(),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5, face = "bold", margin = margin(t = 5), lineheight = 0.75),
               axis.text.y = element_text(size = 9, angle = 0),
               axis.title = element_text(size = 10, color = "black"),
               axis.title.x = element_text(hjust = 0.5, margin = margin(t = 10)),
               axis.title.y = element_text(hjust = 0.45),
               axis.line.x = element_line(),
               axis.line.y = element_line(),
               axis.ticks.x = element_blank(),
               axis.text = element_text(color = "black"),
               legend.position = "right",
               legend.text = element_text(size = 9),
               legend.text.align = 0,
               legend.spacing.x = unit(0.2, "cm"),
               legend.justification = c("left", "bottom"),
               legend.key = element_rect(fill = NA),
               legend.key.size = unit(0.3, "cm"),
               legend.box = "vertical",
               legend.title = element_text(size = 12, face = "bold.italic"),
               panel.grid = element_blank(),
               panel.spacing.x = unit(4.5, "cm"),
               panel.spacing.y = unit(0.25, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 22, margin = margin(l = 0, r = 0)),
               strip.text.y = element_text(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Generating plots in a for loop
for (i in groups$taxa) {
  
  # Generating plot
  p <- get(paste0("plot_", i)) %>%
    mutate(taxon = factor(taxon, levels = levels(taxon))) %>%
    ggplot() +
    geom_bar(aes(x = str_wrap(period, width = 4), y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
    scale_fill_manual(name = i, breaks = levels(get(paste0("plot_", i))$taxon), values = color, labels = get(paste0("names_", i))) +
    {if (i != "total_community") {geom_text(data = get(paste0("whole_", i)),
                                            aes(x = str_wrap(period, width = 4), y = 106, label = paste0(trimws(format(round(mean, 1), nsmall = 1), which = "both"), " ±\n",
                                                                                                         trimws(format(round(sd, 1), nsmall = 1), which = "both"), " %"),
                                                                                 lineheight = 1.00),
              family = "Times", fontface = "bold", hjust = 0.5, size = 4.5)}} +
    labs(x = "Time Period", y = "%") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06)), breaks = seq(0, 100, by = 10)) +
    coord_capped_cart(left = "top") +
    theme +
    {if (i == "Alphaproteobacteria" | i == "Bacteroidota") {theme(axis.title.x = element_blank())}} +
    {if (i == "Bacteroidota" | i == "Cyanobacteria") {theme(axis.title.y = element_blank())}} +
    {if (i == "total_community") {theme(axis.title.y = element_text(hjust = 0.475),
                                        legend.title = element_blank(),
                                        axis.text.x = element_text(size = 9))}}

  # Setting name for object created in loop
  assign(paste0("p", sep = "_", i), p)
  
}

# Combining plots and legends
cowplot <- cowplot :: plot_grid(p_Alphaproteobacteria, p_Bacteroidota,
                                p_Gammaproteobacteria, p_Cyanobacteria,
                                align = "hv", nrow = 2, ncol = 2,
                                labels = c("a", "b", "c", "d"), label_fontfamily = "Times", label_fontface = "bold", label_size = 28) +
  cowplot :: draw_line(x = c(0.375, 0.375), y = c(0.590, 0.627), size = 0.4) +
  cowplot :: draw_label("SAR11 Clade", x = 0.38, y = 0.609, hjust = 0,  fontfamily = "Times", size = 9)

# Saving
ggsave("results/figures/community_bar_plot_month.jpg", p_total_community, width = 210 / 2, height = 297 / 2, units="mm", bg = "white")
ggsave("results/figures/community_bar_plot_month_taxa.jpg", cowplot, width = 210 * 1.25, height = 297, units="mm", bg = "white")
