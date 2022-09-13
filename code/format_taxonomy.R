map.in <- read.table("data/references/tax_slv_ssu_138.txt", header=F,sep="\t", stringsAsFactors=F)
map.in <- map.in[, c(1,3)]
colnames(map.in) <- c("taxlabel", "taxlevel")
taxlevels <- c("root", "domain", "major_clade", "superkingdom", "kingdom", "subkingdom", "infrakingdom", "superphylum", "phylum", "subphylum", "infraphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily", "genus")
taxabb <- c("ro", "do", "mc", "pk", "ki", "bk", "ik", "pp", "ph", "bp", "ip", "pc", "cl", "bc", "ic", "po", "or", "bo", "pf", "fa", "bf", "ge")
tax.mat <- matrix(data="", nrow=nrow(map.in), ncol=length(taxlevels))
tax.mat[,1] <- "root"
colnames(tax.mat) <- taxlevels

outlevels <- c("domain", "phylum", "class", "order", "family", "genus")

for(i in 1:nrow(map.in)) {
  taxname <- unlist(strsplit(as.character(map.in[i,1]), split=';'))
  # Print(taxname);
  
  while ( length(taxname) > 0) {
    # Regex to look for exact match
    
    tax.exp <- paste(paste(taxname,collapse=";"), ";", sep="")
    tax.match <- match(tax.exp,map.in$taxlabel)
    tax.mat[i, map.in[tax.match,2]] <- tail(taxname, 1)
    taxname <- head(taxname, -1)
  }
}

for(i in 1:nrow(tax.mat)) {
  # This fills in the empty gaps by using the closest higher taxonomic level appended with an abbreviation for the current taxonomic level
  # if you don't want this behavior, cut it out
  for(j in 1:ncol(tax.mat)) {
    if(tax.mat[i, j] < 0) { tax.mat[i, j] <- paste(tmptax, taxabb[j], sep="_")}
    else { tmptax <- tax.mat[i, j]}
  }
  
  # This maps the new name to the input taxonomic levels
  map.in[i, "taxout"] <- paste(paste(tax.mat[i, outlevels], collapse=";"), ";", sep="")
}

# Replace spaces with underscores
map.in$taxout <- gsub(" ", "_", map.in$taxout)

# Bring in the old taxonomic levels from SILVA and remap them using the new levels
tax.in <- read.table("data/references/silva.nr_v138.full", header=F, stringsAsFactors=F, sep="\t")
colnames(tax.in) <- c("taxid", "taxlabel")

# Problem corrections observed in the SILVA ARB file (silva.nr_v138.full)
# Following line corrects the Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Mesorhizobium;Mesorhizobium; problem
tax.in$taxlabel <- gsub("Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Mesorhizobium;Mesorhizobium;", "Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Mesorhizobium;", tax.in$taxlabel)

tax.in$id <- 1:nrow(tax.in)

tax.write <- merge(tax.in, map.in, all.x=T, sort=F)
tax.write <- tax.write[order(tax.write$id), ]

# We want to see whether everything has 6 taxonomic level (kingdom to genus)
getDepth <- function(taxonString) {
  initial <- nchar(taxonString)
  removed <- nchar(gsub(";", "", taxonString))
  return(initial-removed)
}

depth <- getDepth(tax.write$taxout)
summary(depth) # Should all be 6 and there should be no NAs
bacteria <- grepl("Bacteria;", tax.write$taxout)
archaea <- grepl("Archaea;", tax.write$taxout)
eukarya <- grepl("Eukaryota;", tax.write$taxout)

tax.write[depth > 6 & bacteria, ] # Good to go
tax.write[depth > 6 & archaea, ]  # Good to go
tax.write[depth > 6 & eukarya, ]  # Good to go

write.table(tax.write[, c("taxid", "taxout")], file="data/references/silva.full_v138.tax", sep="\t", row.names=F, quote=F, col.names=F)
