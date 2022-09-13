############################################################
# Libraries used
############################################################

library(stats)
library(knitr)
library(rmarkdown)
library(tinytex)
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(lemon)

############################################################
# Options for knitr and Rmarkdown rendering
############################################################

## Set working directory for R code chunks (https://bookdown.org/yihui/rmarkdown-cookbook/working-directory.html)
if (require("knitr")) {
    opts_knit$set(root.dir=normalizePath(getwd()))
}

# Avoid false positive error when using knitrinclude_graphics() (knitr release 1.28, https://github.com/yihui/knitr/release/tag/v1.28)
include_graphics = function(...) {
    knitr::include_graphics(..., error = FALSE)
}

############################################################
# Permanently setting the CRAN repository
# Credit: https://www.r-bloggers.com/permanently-setting-the-cran-repository/
############################################################

local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.wu.ac.at/"
  options(repos = r)
})

############################################################
# New key drawing function for legend (ggplot2 issue; https://github.com/tidyverse/ggplot2/issues/2844)
############################################################

# Function that draws key without gap
draw_key_polygon2 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.8, "npc"),
    height = grid::unit(0.8, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# Register new key drawing function, effect is global and persistent
# throughout R session!
GeomBar$draw_key = draw_key_polygon2

############################################################
# Option to keep the auxiliary TeX files when rendering a rmarkdown document
############################################################
options(tinytex.clean = FALSE)
