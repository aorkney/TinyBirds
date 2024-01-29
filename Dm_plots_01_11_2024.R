# The purpose of this script is to produce a combination of plots
# that demonstrate how integration within the head, wing and
# between the wing and trunk, change with increasing body mass.
# This will be achieved by computing major axes of covariance of 
# skeletal proportions within a phylogenetic framework, 
# before determining whether heteroskedasticity (scatter) of the individual
# taxa depends significantly upon body mass. 
# This analysis corresponds to Figure 4 in Orkney & Hedrick 2024. 
#
# This code was written by Dr. Andrew Orkney, parsed by Prof. Brandon Hedrick, 
# and draws upon existing functions from dependent R packages and online contributions from Stackexchange.

library( geomorph ) # 4.0.5
library( ape ) # 5.7-1
library( nlme ) # 3.1-162
library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library( phytools ) # 1.5-1
# Load required packages

setwd('C:/Users/Lab/Documents/Navalon_birds/File S1.Bird_landmarks_Navalon_2022/Clean')
# Set work directory appropriately (you will need to adjust this if you
# are replicating the work undertaken here)

# Load centroid sizes for study birds. 
load('Csize.12.4.2023.RData')
# In our study, we defined this object as an array called 'GPA.Csize'
# This object is derived from 228 birds from across crown Aves. 
# The landmark constellations from which these data were extracted 
# are available at: https://doi.org/10.1038/s41586-022-05372-y
# (Navalon et al., 2022)

# Load phylogenetic tree of birds
load('tree.12.4.2023.RData')
# Pruned tree produced as described in Navalon et al., 2022
# It is a merger of Prum and Oliveros
# The object is a phylogeny called 'pruned.tree'.

# Load bird masses
load('masses.12.4.2023.RData')
# The object is a named numeric vector called 'masses'.
# The names will not yet match those on the phylogeny.

den.tree <- as.dendrogram(force.ultrametric(pruned.tree))
seg.tree<-dendro_data(den.tree)$segments
# Prepare the phylogeny for phylogram plot production
# The phylogeny differed slightly from an ultrametric topology.

# https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2
# I have borrowed code that was discussed in this online forum, in order to construct a plot of
# the avian phylogeny.

library(zoo) # 1.8-12
library(dplyr) # 1.1.1
# These packages will be required to plot the phylogeny.

cut <- 7
# We will colour 7 sub-clades of the avian phylogeny, to help orientate the reader.
dendr <- dendro_data(as.dendrogram(force.ultrametric(pruned.tree))) 
clust <- cutree(force.ultrametric(pruned.tree), k= cut)
clust.df <- data.frame(label = names(clust), cluster = clust)
# A dataframe with the constituent clustering information has been produced.



height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
cut.height <- mean(c(height[cut], height[cut-1]))
dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
   dendr$segments$y > cut.height, 1, 2)
dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)
# These lines divide the dendrogram into the roots and the branches that comprise the 20 sub-clades.

dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
change <- which(dendr$segments$cluster == 1)
for (i in 1:cut) dendr$segments$cluster[change[i]] = i + 1
dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1, 
             ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
dendr$segments$cluster <- na.locf(dendr$segments$cluster) 
# These lines ascribe numbers to the sub-clades

clust.df$label <- factor(clust.df$label, levels = levels(dendr$labels$label))
clust.df <- arrange(clust.df, label)
clust.df$cluster <- factor((clust.df$cluster), levels = unique(clust.df$cluster), labels = (1:cut) + 1)
dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
# These lines assure consistent numbering between segment$cluster and label$cluster

createAngleHJustCols <- function(labeldf) {        
    nn <- length(labeldf$y)
    halfn <- floor(nn/2)
    firsthalf <- rev(90 + seq(0,360, length.out = nn))
    secondhalf <- rev(-90 + seq(0,360, length.out = nn))
    angle <- numeric(nn)
    angle[1:halfn] <- firsthalf[1:halfn]
    angle[(halfn+1):nn] <- secondhalf[(halfn+1):nn]

    hjust <- numeric(nn)
    hjust[1:halfn] <- 0
    hjust[(halfn+1):nn] <- 1

    return(list(angle = angle, hjust = hjust))
}
# This is a function that I borrowed from Stackexchange: https://stackoverflow.com/questions/38034663/rotate-labels-for-ggplot-dendrogram
# The function rotates labels on fan phylogenies. 

lab.dat<-dendro_data(den.tree)$labels
lab.dat$label<- gsub('_.*','',lab.dat$label)
# Extract text to label genera 


# This subplot should be decorated with labelled subclades
# I am going to achieve this by drawing arcs around the periphery of the phyloplot, 
# labelled with major avian subgroups 

setwd('C:/Users/Lab/Documents/Navalon_birds/File S1.Bird_landmarks_Navalon_2022')
# The user will need to change this directory accordinly 
metadata <- read.csv('Navalon_metadata.csv')
# Reading metadata from Navalon et al., 2022 (available at previously specified DOI)
# 38 Orders are represented, so a coarser division is 
# probably necessary
# Let's use the major avian clades in Prum et al., 2015 

Accipitrimorphae <- range(grep( 'Accipitriformes', metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]  ))

ToMatch <- c('Aegotheliformes', 'Caprimulgiformes', 'Apodiformes', 'Podargiformes', 'Steatornithiformes')
Strisores <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))])) 

ToMatch <- c('Anseriformes', 'Galliformes') 
Galloanserae <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))

ToMatch <-c('Bucerotiformes', 'Coliiformes', 'Coraciiformes', 'Leptosomiformes', 'Piciformes', 'Trogoniformes')
Coraciimorphae <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))

ToMatch <-'Gruiformes'
Gruiformes <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))

ToMatch <- c('Cariamiformes', 'Falconiformes') 
Australaves <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))
Australaves[2]<-Australaves[2]+6

ToMatch <-c('Charadriiformes', 'Ciconiiformes', 'Eurypygiformes', 'Gaviiformes', 'Pelecaniformes', 'Phaethontiformes', 'Phoenicopteriformes', 'Podicipediformes', 'Procellariiformes', 'Sphenisciformes', 'Suliformes')
Aequorlithornithes <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))

ToMatch <-c('Columbiformes', 'Cuculiformes', 'Mesitornithiformes', 'Otidiformes', 'Pterocliformes')
Columbaves <- range(grep(paste(ToMatch,collapse="|"), metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]))

# Opisthocomiformes, Psittaciformes --> Inopinaves (near Australavian grade)

Opisthocomiformes <- grep('Opisthocomiformes',metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))])

levels(clust.df$cluster)<-c(levels(clust.df$cluster),'9')
clust.df$cluster[Opisthocomiformes]<-'9'

Passeriformes <- range(grep( 'Passeriformes', metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]  )[-c(1:4)])
# I made a snip here because the first entry 'Sarothrura' is a Gruiform, and not a Passeriform 
# This is a mistake in Navalon et al., 2022's metadata

# Strigiformes --> Strigiformes (near Coraciimorphae)

Strigiformes <- grep('Strigiformes',metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))])

Palaeognathae <- range(grep( 'Tinamiformes', metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]  ))






library(geomtextpath) # version 0.1.1
# This package will be used to curve the text around the phylogram

Path_Palaeognathae <- data.frame( x=Palaeognathae+c(-0,6), y=c(-2.5,-2.5))
Path_Galloanserae <- data.frame(x=Galloanserae, y=c(-1.5,-1.5))
Path_Strisores <- data.frame(x=Strisores, y=c(-1.5,-1.5))
Path_Columbaves <- data.frame(x=Columbaves, y=c(-1.5,-1.5))
Path_Gruiformes <- data.frame(x=Gruiformes, y=c(-1.5,-1.5))
Path_Aequorlithornithes <- data.frame(x=Aequorlithornithes, y=c(-1.5,-1.5))
Path_Accipitrimorphae <- data.frame(x=Accipitrimorphae, y=c(-1.5,-1.5))
Path_Coraciimorphae <- data.frame(x=Coraciimorphae, y=c(-1.5, -1.5))
Path_Passeriformes <- data.frame(x=Passeriformes, y=c(-1.5, -1.5))
Path_near_Passerines <- data.frame(x=Australaves, y=c(-1.5, -1.5))
Path_Opisthocomiformes <- data.frame(x=Opisthocomiformes, y=c(-1.5, -1.5))
Path_Strigiformes <- data.frame(x=Strigiformes, y=c(-1.5, -1.5))
# Curved text paths

# We also desire appealing silhouettes of characteristic taxa 
# These silhouettes are not provided with the code deposition on GitHub
# They cab be downloaded from https://www.phylopic.org/
# As these elements are purely decorative, users seeking to replicate our analyses may run the code without them. 
setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_1_09_2024')
# The user will need to change as appropriate
library(png) # version 0.1-8
# png image processing
library(cowplot)
# grob compilations
Palaeognathae_S <- readPNG("Crypturellus_tataupa.png")
Palaeognathae_S <- grid::rasterGrob(Palaeognathae_S, interpolate=TRUE)
Galloanserae_S <- readPNG("Crax_mitu.png")
Galloanserae_S <- grid::rasterGrob(Galloanserae_S, interpolate=TRUE)
Strisores_S <- readPNG("Archilocus_colubris.png")
Strisores_S <- grid::rasterGrob(Strisores_S, interpolate=TRUE)
Columbaves_S <- readPNG("Columba_palumbus.png")
Columbaves_S <- grid::rasterGrob(Columbaves_S, interpolate=TRUE)
Gruiformes_S <- readPNG("Grus_canadensis.png")
Gruiformes_S <- grid::rasterGrob(Gruiformes_S, interpolate=TRUE)
Aequorlithornithes1_S <- readPNG("Jacana_jacana.png")
Aequorlithornithes1_S <- grid::rasterGrob(Aequorlithornithes1_S, interpolate=TRUE)
Aequorlithornithes2_S <- readPNG("Spheniscus_humboldti.png")
Aequorlithornithes2_S <- grid::rasterGrob(Aequorlithornithes2_S, interpolate=TRUE)
Aequorlithornithes3_S <- readPNG("Sula_dactylatra.png")
Aequorlithornithes3_S <- grid::rasterGrob(Aequorlithornithes3_S, interpolate=TRUE)
Accipitrimorphae_S <- readPNG("Vultur_gryphus.png")
Accipitrimorphae_S <- grid::rasterGrob(Accipitrimorphae_S, interpolate=TRUE)
Coraciimorphae1_S <- readPNG("Upupa_epops.png")
Coraciimorphae1_S <- grid::rasterGrob(Coraciimorphae1_S, interpolate=TRUE)
Coraciimorphae2_S <- readPNG("Alcedo_atthis.png")
Coraciimorphae2_S <- grid::rasterGrob(Coraciimorphae2_S, interpolate=TRUE)
Passeriformes1_S <- readPNG("Acanthisitta_chloris.png")
Passeriformes1_S <- grid::rasterGrob(Passeriformes1_S, interpolate=TRUE)
Passeriformes2_S <- readPNG("Acanthorhynchus_tenuirostris.png")
Passeriformes2_S <- grid::rasterGrob(Passeriformes2_S, interpolate=TRUE)
Passeriformes3_S <- readPNG("Hirundo_rustica.png")
Passeriformes3_S <- grid::rasterGrob(Passeriformes3_S, interpolate=TRUE)
Passeriformes4_S <- readPNG("Turdus_ pilaris.png")
Passeriformes4_S <- grid::rasterGrob(Passeriformes4_S, interpolate=TRUE)


scale<-0.1
# Visual parameter for phyogeny scaling in plots
size <- 4
# Visual parameter for element size scaling in plots 


dendr$segments$cluster[which(dendr$segments$cluster==8)] <- '1' 
dendr$segments$cluster[390] <- '8' 
# The Hoazin is identified as line 390

clust[Opisthocomiformes]<- 8

dendr$segments$cluster[min(4*(Strigiformes-0.5)):max(4*(Strigiformes+0.5))] <- '9' 
dendr$segments$cluster[(min(4*(Strigiformes-0.5)):max(4*(Strigiformes+0.5)))[which(dendr$segments$yend[min(4*(Strigiformes-0.5)):max(4*(Strigiformes+0.5))] >118)]] <- '1' 
clust[min(Strigiformes):(max(Strigiformes))] <- 9

dendr$segments$cluster[min(4*(Accipitrimorphae-0.5)):max(4*(Accipitrimorphae+0.5))] <- '10' 
dendr$segments$cluster[(min(4*(Accipitrimorphae-0.5)):max(4*(Accipitrimorphae+0.5)))[which(dendr$segments$yend[min(4*(Accipitrimorphae-0.5)):max(4*(Accipitrimorphae+0.5))] >120)]] <- '1' 
clust[min(Accipitrimorphae):(max(Accipitrimorphae))] <-10

dendr$segments$cluster[min(4*(Coraciimorphae-0.5)):max(4*(Coraciimorphae-0.5))] <- '11' 
dendr$segments$cluster[(min(4*(Coraciimorphae-0.5)):max(4*(Coraciimorphae-0.5)))[which(dendr$segments$yend[min(4*(Coraciimorphae-0.5)):max(4*(Coraciimorphae-0.5))] >122)]] <- '1' 
clust[min(Coraciimorphae):(max(Coraciimorphae))]<-11

# Before we set passeriformes, 
# we need to set 'Australaves', so that we can get the 
dendr$segments$cluster[min(4*(Australaves-1)):length(dendr$segments$cluster)] <- '12' 
clust[min(Australaves):(max(Australaves))] <- 12
dendr$segments$cluster[min(4*(Passeriformes-1)):length(dendr$segments$cluster)] <- '13' 
clust[min(Passeriformes):(max(Passeriformes))] <- 13

dendr$segments$line[which(dendr$segments$cluster==1)]<-1
#dendr$segments$cluster[c(422,423,425,429,430)]<-'1'
# The line thicknesses of these lines also need to be set to 1

#dendr$segments$cluster[c(395,397,403,407,410)]<-'1'
# The line thicknesses of these lines also need to be set to 1

# Definining additional breaks within the phylogeny, which are to be coloured in the final plot.
# These are complex to define because some of them represent paraphyletic assmeblages (e.g. near Passerines)

cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AA6997", "#885987", "#664977", "#443967",  "#222957")
# A colour-blind friendly palette will be used to colour the major avian groups in the phylogeny

phylogram<-
ggplot()+
	geom_textpath(data=Path_Palaeognathae,size  = size, label = 'Palaeognathae', aes(x=c(0,12),y=y-0.25), linecolour='white')+
	geom_segment(linewidth=2, aes(x=2, xend=2,y=-0.5, yend=-2.0), colour='black')+
	geom_segment(data= Path_Palaeognathae,aes(x=x[1], xend=x[2]-6, y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Galloanserae,size  = size, label = 'Galloanserae', aes(x=x,y=y), linecolour='white')+
	geom_segment(data= Path_Galloanserae,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Strisores,size  = size, label = 'Strisores', aes(x=x,y=y), linecolour='white')+
	geom_segment(data= Path_Strisores,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Columbaves,size  = size, label = 'Columbaves', aes(x=x,y=y), linecolour='white')+
	geom_segment(data= Path_Columbaves,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Gruiformes,size  = size, label = 'Gruiformes', aes(x=c(47,55),y=y), linecolour='white')+
	geom_segment(data= Path_Gruiformes,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Aequorlithornithes,size  = size, label = 'Aequorlithornithes', aes(x=c(56,90),y=y), linecolour='white')+
	geom_segment(data= Path_Aequorlithornithes,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Accipitrimorphae,size  = size, label = 'Accipitrimorphae', aes(x=c(93,103),y=y-1), linecolour='white')+
	geom_segment(linewidth=2, aes(x=102, xend=102,y=-0.5, yend=-2.0), colour='black')+
	geom_segment(data= Path_Accipitrimorphae,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Coraciimorphae,size  = size, label = 'Coraciimorphae', aes(x=c(113,133),y=y), linecolour='white')+
	geom_segment(data= Path_Coraciimorphae,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Passeriformes,size  = size, label = 'Passeriformes', aes(x=x,y=y), linecolour='white')+
	geom_segment(data= Path_Passeriformes,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_near_Passerines,size  = size, label = 'near Passerines', aes(x=x,y=y), linecolour='white')+
	geom_segment(data= Path_near_Passerines,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Opisthocomiformes ,size  = size, label = 'Opisthocomiformes', aes(x=c(86,96),y=y), linecolour='white')+
	geom_segment(data= Path_Opisthocomiformes ,aes(x=97.7, xend=98.5, y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_textpath(data=Path_Strigiformes ,size  = size, label = 'Strigiformes', aes(x=c(105,113),y=y), linecolour='white')+
	geom_segment(data= Path_Strigiformes ,aes(x=x[1], xend=x[2], y=-0.5, yend=-0.5), linewidth=2)+ 
	geom_segment(
	data=segment(dendr),
  	aes(x = x, y = y *scale, xend = xend, yend = yend *scale, size=factor(line), colour=factor(cluster)), 
      lineend = "square", show.legend = FALSE) + 
  	#scale_colour_manual(values = c("grey60", rev(hcl.colors(13,palett='Hawaii')))) + # Alternative colour scale
	scale_colour_manual(values = c("grey60", cbp[c(9,11,12,1,7,2,5,4,6,3,8,10)] )) +
	scale_size_manual(values = c(.1, 1.1)) +
	#geom_text(data= lab.dat, aes(x=x, y=y, label=gsub('_.*','',label)), angle=(createAngleHJustCols(lab.dat)[['angle']]), hjust=(createAngleHJustCols(lab.dat)[['hjust']]),size = 6 / .pt )+
	# Individual taxon labels may be plotted by Users who seek to replicate and inspect our analyses. 
	scale_y_reverse(expand = c(0.2, 0)) + 
	coord_polar()+
	expand_limits(x=0)+
	theme(legend.position='bottom',
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
      )


h <- ggdraw(phylogram)
h<- h + draw_grob(Palaeognathae_S, 0.48, 0.85, 0.08, 0.08)+
 draw_grob(Galloanserae_S, 0.62, 0.815, 0.08, 0.08)+
 draw_grob(Strisores_S, 0.72, 0.75, 0.08, 0.08)+
 draw_grob(Columbaves_S, 0.79, 0.65, 0.08, 0.08)+
 draw_grob(Gruiformes_S, 0.835, 0.52, 0.11, 0.11)+
 draw_grob(Aequorlithornithes1_S, 0.82, 0.38, 0.11, 0.11)+
 draw_grob(Aequorlithornithes2_S, 0.8, 0.29, 0.1, 0.1)+
 draw_grob(Aequorlithornithes3_S, 0.70, 0.14, 0.11, 0.11)+
 draw_grob(Accipitrimorphae_S, 0.60, 0.07, 0.10, 0.10)+
 draw_grob(Coraciimorphae1_S, 0.47, 0.07, 0.1, 0.1)+
 draw_grob(Coraciimorphae2_S, 0.35, 0.10, 0.08, 0.08)+
 draw_grob(Passeriformes1_S, 0.15, 0.22, 0.09, 0.09)+
 draw_grob(Passeriformes2_S, 0.08, 0.5, 0.11, 0.11)+
 draw_grob(Passeriformes3_S, 0.17, 0.72, 0.11, 0.11)+
 draw_grob(Passeriformes4_S, 0.31, 0.83, 0.08, 0.08)

# Now that we have the silhouettes, colours and taxon labels we want to be able to produce subplots in R as described at the start of the code.
# User seeking to replicate our analysis without phylograms can simply run the line 'h<-ggdraw(phylogram)' and ignore the subsequent 'draw_grob' commands. 

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute landmark constellation centroid sizes, corrected for allometric scaling.
# We want to account for 'allometry' (size-dependent scaling between species), because it could itself be a major driver of trait integration.
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # Input arguments of shape-data array, mass vector, phylogeny, taxa and bones of interest
	allometry.Csize <- list() # Dumby variable defined to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to taxa of interest
	for(i in 1:length(array[bones]) ){ # For each bone of interest 
				lambda <- phylosig(tree= newphy, x=log10( array[[ i ]][taxa]),  method='lambda')[[1]] # Initial estimate of phylogenetic autocorrelation
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Organise the data into a dataframe
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( lambda, phy=newphy, form= ~ species  ), data=df ) # Compute allometric model
	}
	names(allometry.Csize) <- names(array[bones]) # Check correspondance of allometric models with bone names that describe them 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract allometric model residuals 
	return(residual_Csize) # Return residuals 
}

# Define a function to compute the major axis of trait covariance between two sets of skeletal proportions. 
m.axis.resid <- function( array, masses, bone1, bone2, phylogeny, taxa ){ # This function takes shape data, a mass vector, bone names, phylogeny and taxa of interest inputs
	bones<-c(bone1,bone2) # Define a pairwise combination of bones 
	r.Csize <- get.residual.Csize( array = array, masses = masses, phylogeny = phylogeny, taxa= taxa,bones= bones ) # Compute residual centroid size
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune phylogeny to taxa of interest 
	two.block <- phylo.integration(r.Csize[[bone1]],r.Csize[[bone2]],phy=newphy) # Compute 2 block partial least squares (integration) 
	residuals<- prcomp(cbind(two.block$XScores[,1],two.block$YScores[,1]))$x[,2] # Find residuals from the major axis of covaraince 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals
	return(list) # Return residuals 
}


sample.bins <- 10^seq(log10(min(masses)),log10(max(masses)),length.out=10)
# Define breaks across the vector of avian body masses to sample within
# This is necessary to achieve a Gaussian-like distribution of sampled masses.
# The total distribution has a right skewness. 


bone1<-'carpometacarpus'
bone2<-'humerus'
# Let us consider whether integration between the carpometacarpus and humerus changes with body mass across birds 

k<-1
p.list <- list()
while(length(p.list) < 100){

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses >= sample.bins[i] & masses <= sample.bins[i+1]), 10,replace=T))) # Gaussian-like subsample
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} ) # Find major axis 
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses[sample.list]) )
		p.list[[k]] <- car::ncvTest(model, ~log10(masses[sample.list]))$p #3.1.1 # Perform Breusch-Pagan test for unequal scatter of residuals 
		k <- k+1
	}
	print(k)
}

# Perform above analysis again, but in a scrambled dataset that ipso-facto exhibits no patterns:
k<-1
p.list.permuted <- list()
while(length(p.list.permuted) < 100){

	birds.permuted <- GPA.Csize
	
	for(i in 1:length(birds.permuted)){
		birds.permuted[[i]]<- sample(birds.permuted[[i]])
		names(birds.permuted[[i]]) <- names(masses)
	}
	
	masses.permuted <- masses
	names(masses.permuted)<-sample(names(masses))

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses.permuted >= sample.bins[i] & masses.permuted <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( birds.permuted, masses.permuted, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses.permuted[sample.list]) )
		p.list.permuted[[k]] <- car::ncvTest(model, ~log10(masses.permuted[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}

dv <- data.frame(cbind(rep(c('real','permuted'),each=100), c(-log10(unlist(p.list)),-log10(unlist(p.list.permuted)))))
dv$X2<-as.numeric(dv$X2)
# Compile results into a dataframe 



wing.inset<- # 'Flame plot' which illustrates the distribution of p-values for an unequal scatter test (difference in intergation over body mass) compared to a null
# distribution in grey. 
ggplot()+
geom_violin( data=dv[101:200,], aes(x=1, y=X2), fill='darkgrey',colour='NA')+
geom_violin( data=dv[1:100,], aes(x=1, y=X2), fill=NA,colour='black',lwd=1)+
#geom_sina( data=dv[1:100,], aes(x=1, y=X2), alpha=0.5,colour='black', shape=21,size=2)+
labs(y=expression(paste(-log[10],'(',italic(p),')',' ')) )+
theme(axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.y=element_text(size=15),
				axis.ticks.x=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
      			panel.background=element_blank(),
      			panel.border=element_blank(),
      			panel.grid.major=element_blank(),
      			panel.grid.minor=element_blank(),
      			plot.background=element_blank(),
				)

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )
# Find the residuals from the major axis relating bones 1 and 2.
# We will need to visualise the relationship to determine its sense. 

library(caTools) # 1.18.2
# Package with a running quantile function
# We will now begin to visualise the relationship.
# We need to load this package to use running quantile functions for plotting sigma contour intervals. 

# Compute 1 and 2 sigma confidence intervals for the data distribution
# A running window of k=30 is applied to comptue quantiles 
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

library(ggrepel) # 0.9.3
# This package can be applied to minimise label overlap.

 axis.title <- 15
 axis.text <- 10
 point.size <- 3
 line.width <- 1/2
 genus.text<-2.5
# Various plotting parameters

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
# Compile a dataframe for plotting

outside <-  rev(order(abs(model$residuals)))[1:9] 
# Identify outliers 

wing.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1/2, data=df, col='black', pch=21, show.legend=F)+
  	#scale_fill_manual(values = c( rev(hcl.colors(13,palett='Hawaii')))) +
scale_fill_manual(values = c( cbp[c(7,2,5,4,6,3,8,10,9,11,12,1)] )) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 11 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='humerus-carpometacarpus')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))
wing.plot.inset <- wing.plot+
draw_plot(wing.inset, x = max(df$V1)*0.77, y = min(df$V2), width = max(df$V1)*0.3, height = diff(range(df$V2))/2 )
# This plot shows that the carpometacarpus and humerus tend to evolve more indpendently of one another in the smallest birds,
# whereas they are tightly integrated at a macroevolutionary scale in large birds. 
# The inset flame plot demonstrates that the result is highly significant compared to a null distribution. 

# Hereafter analyses are repeated for other combinations of bones: 

bone1<-'skull'
bone2<-'mandible'

k<-1
p.list <- list()
while(length(p.list) < 100){

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses >= sample.bins[i] & masses <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses[sample.list]) )
		p.list[[k]] <- car::ncvTest(model, ~log10(masses[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}


k<-1
p.list.permuted <- list()
while(length(p.list.permuted) < 100){

	birds.permuted <- GPA.Csize
	
	for(i in 1:length(birds.permuted)){
		birds.permuted[[i]]<- sample(birds.permuted[[i]])
		names(birds.permuted[[i]]) <- names(masses)
	}
	
	masses.permuted <- masses
	names(masses.permuted)<-sample(names(masses))

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses.permuted >= sample.bins[i] & masses.permuted <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( birds.permuted, masses.permuted, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses.permuted[sample.list]) )
		p.list.permuted[[k]] <- car::ncvTest(model, ~log10(masses.permuted[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}

dv <- data.frame(cbind(rep(c('real','permuted'),each=100), c(-log10(unlist(p.list)),-log10(unlist(p.list.permuted)))))
dv$X2<-as.numeric(dv$X2)



head.inset<-
ggplot()+
geom_violin( data=dv[101:200,], aes(x=1, y=X2), fill='darkgrey',colour='NA')+
geom_violin( data=dv[1:100,], aes(x=1, y=X2), fill=NA,colour='black',lwd=1)+
#geom_sina( data=dv[1:100,], aes(x=1, y=X2), alpha=0.5,colour='black', shape=21,size=2)+
labs(y=expression(paste(-log[10],'(',italic(p),')',' ')) )+
theme(axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.y=element_text(size=15),
				axis.ticks.x=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
      			panel.background=element_blank(),
      			panel.border=element_blank(),
      			panel.grid.major=element_blank(),
      			panel.grid.minor=element_blank(),
      			plot.background=element_blank(),
				)

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )

# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 


df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))

outside <-  rev(order(abs(model$residuals)))[1:6] 

head.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1/2, data=df, col='black', pch=21, show.legend=F)+
  	#scale_fill_manual(values = c( rev(hcl.colors(13,palett='Hawaii')))) +
scale_fill_manual(values = c( cbp[c(7,2,5,4,6,3,8,10,9,11,12,1)] )) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 11 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='cranium-mandible')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))
head.plot.inset <- head.plot+
draw_plot(head.inset, x = max(df$V1)*0.77, y = min(df$V2), width = max(df$V1)*0.3, height = diff(range(df$V2))/2 )

# This plot shows that there is little visual evidence of any change in integration between the cranium and mandible as a function of body mass
# across birds. There is a slightly pinching around clade passeriformes, because of high sampling within this clade (vis a vis central limits theorem). 
# The inset flameplot shows that the real and null distributions exhibit substantial overlap, showing that this dataset is not substantially more structured
# than a random dataset. 

bone2<-'scapula'
bone1<-'sternum'

k<-1
p.list <- list()
while(length(p.list) < 100){

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses >= sample.bins[i] & masses <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses[sample.list]) )
		p.list[[k]] <- car::ncvTest(model, ~log10(masses[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}


k<-1
p.list.permuted <- list()
while(length(p.list.permuted) < 100){

	birds.permuted <- GPA.Csize
	
	for(i in 1:length(birds.permuted)){
		birds.permuted[[i]]<- sample(birds.permuted[[i]])
		names(birds.permuted[[i]]) <- names(masses)
	}
	
	masses.permuted <- masses
	names(masses.permuted)<-sample(names(masses))

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses.permuted >= sample.bins[i] & masses.permuted <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( birds.permuted, masses.permuted, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses.permuted[sample.list]) )
		p.list.permuted[[k]] <- car::ncvTest(model, ~log10(masses.permuted[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}

dv <- data.frame(cbind(rep(c('real','permuted'),each=100), c(-log10(unlist(p.list)),-log10(unlist(p.list.permuted)))))
dv$X2<-as.numeric(dv$X2)



trunk.inset<-
ggplot()+
geom_violin( data=dv[101:200,], aes(x=1, y=X2), fill='darkgrey',colour='NA')+
geom_violin( data=dv[1:100,], aes(x=1, y=X2), fill=NA,colour='black',lwd=1)+
#geom_sina( data=dv[1:100,], aes(x=1, y=X2), alpha=0.5,colour='black', shape=21,size=2)+
labs(y=expression(paste(-log[10],'(',italic(p),')',' ')) )+
theme(axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.y=element_text(size=15),
				axis.ticks.x=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
      			panel.background=element_blank(),
      			panel.border=element_blank(),
      			panel.grid.major=element_blank(),
      			panel.grid.minor=element_blank(),
      			plot.background=element_blank(),
				)

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )

# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))


outside <-  rev(order(abs(model$residuals)))[1:7] 


trunk.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1/2, data=df, col='black', pch=21, show.legend=F)+
  	#scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
scale_fill_manual(values = c( cbp[c(7,2,5,4,6,3,8,10,9,11,12,1)] )) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 11 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='scapula-sternum')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))
trunk.plot.inset <- trunk.plot+
draw_plot(trunk.inset, x = max(df$V1)*0.12, y = min(df$V2), width = max(df$V1)*0.3, height = diff(range(df$V2))/2 )
# This plot shows some ambiguous evidence that suggests it is possible that the sternum and scapula are more tightly integrated at a macroevolutionary scake
# within small birds, compared to large birds. 


bone1<-'carpometacarpus'
bone2<-'sternum'

k<-1
p.list <- list()
while(length(p.list) < 100){

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses >= sample.bins[i] & masses <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses[sample.list]) )
		p.list[[k]] <- car::ncvTest(model, ~log10(masses[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}


k<-1
p.list.permuted <- list()
while(length(p.list.permuted) < 100){

	birds.permuted <- GPA.Csize
	
	for(i in 1:length(birds.permuted)){
		birds.permuted[[i]]<- sample(birds.permuted[[i]])
		names(birds.permuted[[i]]) <- names(masses)
	}
	
	masses.permuted <- masses
	names(masses.permuted)<-sample(names(masses))

	sample.list <- list()
	for(i in 1:(length(sample.bins)-1) ){
		sample.list[[i]] <- unique(names(sample( which(masses.permuted >= sample.bins[i] & masses.permuted <= sample.bins[i+1]), 10,replace=T)))
	}
	sample.list <-unlist(sample.list)
	hist(log10(masses[sample.list]))
	
	output<-tryCatch( {m.axis.resid( birds.permuted, masses.permuted, bone1, bone2, pruned.tree, sample.list)},
    	error = function(e){output <-1} )
	
	if(length(output)>1 ){
		model <- lm(output[[2]][match(sample.list, names(output[[2]]))] ~log10(masses.permuted[sample.list]) )
		p.list.permuted[[k]] <- car::ncvTest(model, ~log10(masses.permuted[sample.list]))$p #3.1.1
		k <- k+1
	}
	print(k)
}

dv <- data.frame(cbind(rep(c('real','permuted'),each=100), c(-log10(unlist(p.list)),-log10(unlist(p.list.permuted)))))
dv$X2<-as.numeric(dv$X2)



cross.inset<-
ggplot()+
geom_violin( data=dv[101:200,], aes(x=1, y=X2), fill='darkgrey',colour='NA')+
geom_violin( data=dv[1:100,], aes(x=1, y=X2), fill=NA,colour='black',lwd=1)+
#geom_sina( data=dv[1:100,], aes(x=1, y=X2), alpha=0.5,colour='black', shape=21,size=2)+
labs(y=expression(paste(-log[10],'(',italic(p),')',' ')) )+
theme(axis.text.y=element_text(size=15,colour='black',face = "bold"),
				axis.title.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.y=element_text(size=15),
				axis.ticks.x=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
      			panel.background=element_blank(),
      			panel.border=element_blank(),
      			panel.grid.major=element_blank(),
      			panel.grid.minor=element_blank(),
      			plot.background=element_blank(),
				)

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )


# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))

outside <-  rev(order(abs(model$residuals)))[1:10] 




wing.trunk.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1/2, data=df, col='black', pch=21, show.legend=F)+
  	#scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
scale_fill_manual(values = c( cbp[c(7,2,5,4,6,3,8,10,9,11,12,1)] )) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 11 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='carpometacarpus-sternum')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))
cross.plot.inset <- wing.trunk.plot+
draw_plot(cross.inset, x = max(df$V1)*0.06, y = max(df$V2)/3.3, width = max(df$V1)*0.3, height = diff(range(df$V2))/2 )

# This plot provides convincing evidence that the carpometacarpus and sternum are more strongly integrated at a macroevolutionary scale in small birds, 
# and that they evolve more independently in larger birds. 



library(ggpubr) # 0.6.0 
library(grid) # base
# These packages contain functions that we can use to combine subplots:

quadrate <- ggarrange(ncol=2,nrow=2,wing.plot.inset,head.plot.inset,trunk.plot.inset,cross.plot.inset, labels=c('b','c','d','e'),font.label=list(size=25) , hjust=.7, vjust= 3)
# Arange righthand subplots

annotated<- annotate_figure(quadrate, left = textGrob(expression(paste( '',italic('D'[m]),'' )), rot = 90, vjust = 0.5,  gp = gpar(fontsize=20)),
bottom= textGrob(expression(paste(log[10],'(',italic(mass),')',' [g]')), gp = gpar(fontsize=20)))
# Label subplots

dev.new(width=18,height=9,unit='cm')
ggarrange(h,annotated, labels=c('a'),font.label=list(size=25), hjust=-3, vjust= 3 )
# Produce final plot

setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_1_09_2024')
# Set the work directory to the location where you wish to save the combined plots; User must re-specify 
ggsave(filename='Dm_plot_01_11_2024.pdf')
# Save the plots 


