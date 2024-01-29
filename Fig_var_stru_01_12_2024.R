# This script is part of a project that seeks to determine how variation in body mass structures patterns of 
# interspecific evolutionary covariance in skeletal proportions over the avian skeleton. 
# We should expect locomotory structures, such as wings, to become more integrated as mechanical stress increases, 
# which might have implications for the pathways of evolutionary change that are available to avian subclades exploring
# different ranges of body masses.
#
# The purpose of this script is to produce a plot (Figure 1 in manuscript)
# That illustrates how avian body mass is distributed across phylogeny. 
# (subplot a) and perform tests for phylogenetic signal (Blomberg's K and Pagel's Lambda)
#
# In addition to this, raw variance structure for the different
# bone centroid sizes of the avian skeleton need to be plotted 
# as a function of mass.
# (subplot b)
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

library(geomtextpath) # 0.1.1
# This package will be used to curve the text around the phylogram


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

dendr <- dendro_data(as.dendrogram(force.ultrametric(pruned.tree))) 
basic <- as.dendrogram(force.ultrametric(pruned.tree))
# Prepare the phylogeny for phylogram plot production
# The phylogeny differed slightly from an ultrametric topology.

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

den.tree <- as.dendrogram(force.ultrametric(pruned.tree))
lab.dat<-dendro_data(den.tree)$labels
lab.dat$label<- gsub('_.*','',lab.dat$label)
# Extract text to label genera 

df.plot <- as.data.frame( log10(masses[match(pruned.tree$tip,names(masses))]) )
colnames(df.plot)<-'mass'
# Prepare a dataframe of log-transformed bird masses for plotting

fit<-fastAnc(pruned.tree,log10(masses)[pruned.tree$tip],vars=TRUE,CI=TRUE)
# Estimate the ancestral character states for each node
# We will colour each line according to the node value

n.col <- rep(NA, length(dendr$segments[,1]))
n.col[which(dendr$segments$x==dendr$segments$xend)] <- fit$ace[match(pruned.tree$edge[,1],names(fit$ace))]
n.col[which(dendr$segments$x==dendr$segments$xend)-1] <- fit$ace[match(pruned.tree$edge[,1],names(fit$ace))]
dendr.col<- cbind(dendr$segments, n.col  )
colnames(dendr.col)[5]<-'col'
dendr.col <- data.frame(dendr.col)
# Set information to colour tip and branch colours (these will indicate body mass variety over phylogeny)

BlomK <- phylosig(tree=pruned.tree, x= log10(masses)[pruned.tree$tip], method='K',test=T )[[1]]
# Blomberg's K statistic is 1.53
lambda <- phylosig(tree=pruned.tree, x= log10(masses)[pruned.tree$tip], method='lambda',test=T )[[1]]
# Pagel's lambda is 0.92
# A K > 1 indicates that variance tends to be partitioned among clades more strongly than expected 
# under Brownian motion; clades have therefore diverged more than expected 

scale<-.1
# This defines the scale of the root to tip length of the phylogeny.
size<-4
# Plot element size

# This subplot should be decorated labelled subclades
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

Passeriformes <- range(grep( 'Passeriformes', metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]  )[-c(1:4)])
# I made a snip here because the first entry 'Sarothrura' is a Gruiform, and not a Passeriform 
# This is a mistake in Navalon et al., 2022's metadata

# Strigiformes --> Strigiformes (near Coraciimorphae)

Strigiformes <- grep('Strigiformes',metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))])

Palaeognathae <- range(grep( 'Tinamiformes', metadata$Order[match(lab.dat$label,gsub('_.*', '', metadata$Data))]  ))


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
# Define curved paths for clade labels

# The following will produce a labelled phylogram 
plotA <- 
ggplot()+
	geom_textpath(data=Path_Palaeognathae,size  = size, label = 'Palaeognathae', aes(x=c(0,12),y=y-4.25), linecolour='white')+
	geom_segment(linewidth=2, aes(x=2, xend=2,y=-4, yend=-5.5), colour='black')+
	geom_segment(data= Path_Palaeognathae,aes(x=x[1], xend=x[2]-6, y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Galloanserae,size  = size, label = 'Galloanserae', aes(x=x,y=y-4), linecolour='white')+
	geom_segment(data= Path_Galloanserae,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Strisores,size  = size, label = 'Strisores', aes(x=x,y=y-4), linecolour='white')+
	geom_segment(data= Path_Strisores,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Columbaves,size  = size, label = 'Columbaves', aes(x=x,y=y-4), linecolour='white')+
	geom_segment(data= Path_Columbaves,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Gruiformes,size  = size, label = 'Gruiformes', aes(x=c(47,55),y=y-4), linecolour='white')+
	geom_segment(data= Path_Gruiformes,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Aequorlithornithes,size  = size, label = 'Aequorlithornithes', aes(x=c(56,90),y=y-4), linecolour='white')+
	geom_segment(data= Path_Aequorlithornithes,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Accipitrimorphae,size  = size, label = 'Accipitrimorphae', aes(x=c(93,103),y=y-5), linecolour='white')+
	geom_segment(linewidth=2, aes(x=102, xend=102,y=-4, yend=-6.0), colour='black')+
	geom_segment(data= Path_Accipitrimorphae,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Coraciimorphae,size  = size, label = 'Coraciimorphae', aes(x=c(113,133),y=y-4), linecolour='white')+
	geom_segment(data= Path_Coraciimorphae,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Passeriformes,size  = size, label = 'Passeriformes', aes(x=x,y=y-4), linecolour='white')+
	geom_segment(data= Path_Passeriformes,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_near_Passerines,size  = size, label = 'near Passerines', aes(x=x,y=y-4), linecolour='white')+
	geom_segment(data= Path_near_Passerines,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Opisthocomiformes ,size  = size, label = 'Opisthocomiformes', aes(x=c(86,96),y=y-4), linecolour='white')+
	geom_segment(data= Path_Opisthocomiformes ,aes(x=97.7, xend=98.5, y=-4, yend=-4), linewidth=2)+ 
	geom_textpath(data=Path_Strigiformes ,size  = size, label = 'Strigiformes', aes(x=c(105,113),y=y-4), linecolour='white')+
	geom_segment(data= Path_Strigiformes ,aes(x=x[1], xend=x[2], y=-4, yend=-4), linewidth=2)+ 
	geom_col(data = df.plot,aes(x=1:dim(df.plot)[1], y=rep(-3.5,dim(df.plot)[1])), fill='black', linewidth=1, colour='black' )+
	geom_col(data = df.plot,aes(x=1:dim(df.plot)[1], y=rep(-3.3,dim(df.plot)[1]), fill=mass, col=mass), linewidth=1 )+
	scale_fill_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
	scale_color_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
	geom_segment(data=dendr.col, aes(x = x, y = y *scale, xend = xend, yend = yend *scale ),linewidth=2, col='black')+
	geom_segment(data=dendr.col, aes(x = x, y = y *scale, xend = xend, yend = yend *scale, col=col ),linewidth=1)+
	scale_color_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
	#geom_text(data= lab.dat, aes(x=x, y=y, label=label), angle=createAngleHJustCols(lab.dat)[['angle']], hjust=createAngleHJustCols(lab.dat)[['hjust']],size = 5 / .pt )+
	# Taxon labels (currently disabled) User may wish to view individual taxa and their distribution across this tree
	labs(fill=expression(paste(log[10],'(',italic(mass),')',' [g]')), col= expression(paste(log[10],'(',italic(mass),')',' [g]')))+
	scale_y_reverse(expand = c(0.2, 0)) + 
 	coord_polar()+
	#geom_text( aes(x= c(144), y= c(-20.5)  ) , label= as.expression(bquote(K%~~%.(round(BlomK,digits=2))~','~lambda%~~%.(round(lambda,digits=2)))), size=20 /.pt )+
	#geom_text( aes(x= c(134), y= c(-13.5)  ) , label= as.expression(bquote(K%~~%.(round(BlomK,digits=2))~','~lambda%~~%.(round(lambda,digits=2)))), size=20 /.pt )+
	# Code that can super-impose the phylogenetic signal values on the plot, if desired
	theme(legend.position=c(0.3,0.90),
	legend.direction='horizontal',
	legend.title=element_text(size= 20 ),
	legend.text=element_text(size=15),
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
	#plot.margin=margin(rep(-10,4)),
	legend.box.spacing = unit(0, "pt"),
	plot.margin=unit(c(-0.05,-0.10,-0.05,-0.10), "null"),
      )
# This is a coloured phylogeny plot of body mass over the study bird taxa.


# We also desire appealing silhouettes of characteristic taxa 
# These silhouettes are not provided with the code deposition on GitHub
# They cab be downloaded from https://www.phylopic.org/
# As these elements are purely decorative, users seeking to replicate our analyses may run the code without them. 

setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_1_09_2024')
# The user will need to change as appropriate
library(png)
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

h <- ggdraw(plotA)
h<- h + draw_grob(Palaeognathae_S, 0.48, 0.9, 0.08, 0.08)+
 draw_grob(Galloanserae_S, 0.63, 0.855, 0.08, 0.08)+
 draw_grob(Strisores_S, 0.74, 0.78, 0.08, 0.08)+
 draw_grob(Columbaves_S, 0.83, 0.67, 0.08, 0.08)+
 draw_grob(Gruiformes_S, 0.865, 0.52, 0.11, 0.11)+
 draw_grob(Aequorlithornithes1_S, 0.865, 0.38, 0.11, 0.11)+
 draw_grob(Aequorlithornithes2_S, 0.83, 0.26, 0.1, 0.1)+
 draw_grob(Aequorlithornithes3_S, 0.75, 0.12, 0.11, 0.11)+
 draw_grob(Accipitrimorphae_S, 0.60, 0.03, 0.10, 0.10)+
 draw_grob(Coraciimorphae1_S, 0.47, 0.03, 0.1, 0.1)+
 draw_grob(Coraciimorphae2_S, 0.35, 0.05, 0.08, 0.08)+
 draw_grob(Passeriformes1_S, 0.12, 0.22, 0.09, 0.09)+
 draw_grob(Passeriformes2_S, 0.03, 0.5, 0.11, 0.11)+
 draw_grob(Passeriformes3_S, 0.14, 0.74, 0.11, 0.11)+
 draw_grob(Passeriformes4_S, 0.31, 0.83, 0.08, 0.08)
# Add silhouettes to the plot. 
# If the user is not adding silhouettes, simply run the command 'h<-ggdraw(plotA)' and do not add draw_grob commands

# The next task is to assess whether there is 
# heteroskedasticity in the variance structure of the 
# residual centroid sizes which will be the subject of integration tests.

# Define a sundry function to interract with data objects such as GPA.coords
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute the residual skeletal element centroid cizes, after the affects of allometry are removed.
# We want to account for 'allometry' (size-dependent scaling between species), because it could itself be a major driver of trait integration.
get.residual.Csize <- function( array, masses, phylogeny, taxa ){ # The function takes a shape array, mass vector, phylogeny and list of taxa
	species<-taxa
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to desired taxa
	for(i in 1:length(array) ){ # For each skeletal element 
			lambda <- phylosig(tree= newphy, x=log10( array[[ i ]][taxa]),  method='lambda')[[1]] # Initial estimate of phylogenetic autocorrelation
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  ) # Define a data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( lambda, phy=newphy, form= ~ species ), data=df ) # Compute an allometric model
	}
	names(allometry.Csize) <- names(array) # Ensure names of allometric models match the skeletal elements they describe 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract model residuals
	return(residual_Csize) # Return residuals
} # Conclude function
# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 

# Compute residual centroid sizes by applying function:
residual.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses) )
# Which taxa do we wish to investigate? 
# In this example, we will investigate all birds. 

# An array of residual centroid sizes is now available for investigation. 
# 
# We will compute a proxy of variance over a running window of 30 taxa organised by mass, 
# and then plot this against mass.
# The caTools package has a function 'runquantile' that can be used to produce the 
# the interquartile range of a running bin of 30 birds 
# This will need to be plotted as four subplots- one for bones within each module, 
# and each will require their own legend to distinguish line types. 
# They can share a common axis for mass, but they will have different axes
# for the range of residuals 

library(caTools) # 1.18.2

inter.quart <- function(input,k){
	return(runquantile(x=input[order(masses)],k=k,probs=.75,endrule='quantile')-
	runquantile(x=input[order(masses)],k=k,probs=.25,endrule='quantile'))
}
# Compute the interquartile range (range between first and third quartile of variance in skeletal proportions)
# over body mass.
# k = number of values included in each interquartile range
# input = ordered masses from lowest to highest

df<- cbind(rep(masses[order(masses)],length(residual.Csize)), 
rep(names(residual.Csize), each=length(masses)),
unlist(lapply(residual.Csize, inter.quart,k=30))
)
df<-as.data.frame(df)
colnames(df)<-c('mass','bone','range')
df$mass<-as.numeric(df$mass)
df$range<-as.numeric(df$range)
# Compiled results into a dataframe.

# Produce plots of variance structure with body mass for individual modules in the avian body plan.
# Head
bi<-
ggplot( data= df[which(df$bone=='skull' | df$bone=='mandible'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=2 )+
 geom_curve(  aes(xend = 2.1, yend = 0.2, x = 1.6, y = 0.25), arrow = arrow(length = unit(0.03, "npc") ), colour= '#E69F00')+
geom_text(aes(x=1.6,y=0.26),label='mandible', colour= '#E69F00', size=5)+
 geom_curve(  aes(xend = 3.2, yend = 0.1, x = 3.6, y = 0.07), arrow = arrow(length = unit(0.03, "npc") ), colour= '#D55E00')+
geom_text(aes(x=3.6,y=0.06),label='cranium', colour= '#D55E00', size=5)+
scale_colour_manual(values=c('#E69F00','#D55E00'))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	legend.position='none',
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Wing
bii<-
ggplot( data= df[which(df$bone=='humerus' | df$bone=='ulna' | df$bone=='radius' | df$bone=='carpometacarpus'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=2 )+
scale_colour_manual(values=c('#AA9E73','black','#009E73','darkgreen'))+
 geom_curve(  aes(xend = 2.1, yend = 0.1, x = 2.6, y = 0.05), arrow = arrow(length = unit(0.03, "npc") ), colour= 'black')+
geom_text(aes(x=2.6,y=0.04),label='humerus', colour= 'black', size=5)+
 geom_curve(  aes(xend = 3.8, yend = 0.08, x = 3.4, y = 0.04), arrow = arrow(length = unit(0.03, "npc") ), colour= '#AA9E73')+
geom_text(aes(x=3.4,y=0.03),label='carpometacarpus', colour= '#AA9E73', size=5)+
 geom_curve(  aes(xend = 2.3, yend = 0.27, x = 1.5, y = 0.27), arrow = arrow(length = unit(0.03, "npc") ), colour= 'darkgreen')+
geom_text(aes(x=1.5,y=0.28),label='ulna', colour= 'darkgreen', size=5)+
 geom_curve(  aes(xend = 2.55, yend = 0.295, x = 3.1, y = 0.29), arrow = arrow(length = unit(0.03, "npc") ), colour= '#009E73')+
geom_text(aes(x=3.1,y=0.28),label='radius', colour= '#009E73', size=5)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	legend.position='none',
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Trunk
biii<-
ggplot( data= df[which(df$bone=='scapula' | df$bone=='coracoid' | df$bone=='sternum' | df$bone=='synsacrum'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=2 )+
scale_colour_manual(values=c('darkblue','#56B4E9','#0072B2','blue'))+
 geom_curve(  aes(xend = 1.95, yend = 0.14, x = 2.8, y = 0.12), arrow = arrow(length = unit(0.03, "npc") ), colour= 'darkblue')+
geom_text(aes(x=2.8,y=0.115),label='coracoid', colour= 'darkblue', size=5)+
 geom_curve(  aes(xend = 1.5, yend = 0.09, x = 1, y = 0.12), arrow = arrow(length = unit(0.03, "npc") ), colour= 'blue')+
geom_text(aes(x=1,y=0.125),label='synsacrum', colour= 'blue', size=5)+
 geom_curve(  aes(xend = 4, yend = 0.05, x = 3.3, y = 0.03), arrow = arrow(length = unit(0.03, "npc") ), colour= '#56B4E9')+
geom_text(aes(x=3.3,y=0.025),label='scapula', colour= '#56B4E9', size=5)+
 geom_curve(  aes(xend = 3.75, yend = 0.095, x = 3.5, y = 0.12), arrow = arrow(length = unit(0.03, "npc") ), colour= '#0072B2')+
geom_text(aes(x=3.5,y=0.125),label='sternum', colour= '#0072B2', size=5)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	legend.position='none',
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Leg
biv<-
ggplot( data= df[which(df$bone=='femur' | df$bone=='tibiotarsus' | df$bone=='tarsometatarsus'), ] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=2 )+
scale_colour_manual(values=c('darkred','hotpink','#CC79A7'))+
 geom_curve(  aes(xend = 3, yend = 0.35, x = 1.5, y = 0.35), arrow = arrow(length = unit(0.03, "npc") ), colour= 'hotpink')+
geom_text(aes(x=1.5,y=0.365),label='tarsometatarsus', colour= 'hotpink', size=5)+
geom_text(aes(x=3.7,y=0.15),label='tibiotarsus', colour= '#CC79A7', size=5)+
geom_text(aes(x=3,y=0.08),label='femur', colour= 'darkred', size=5)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	legend.position='none',
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

library(ggpubr) # 0.6.0
library(grid) # base
# These packages contain functions that we can use to combine subplots:

quadrate <- ggarrange(bii,bi,biii,biv,ncol=2, nrow=2, labels=c('Wing','Head','Trunk','Leg'), hjust=-1.5)
annotated<- annotate_figure(quadrate, left = textGrob(expression(paste('Interquartile range ',log[10],'(',italic(rCsize),')',' ')), rot = 90, vjust = 0.5,  gp = gpar(cex = 1.8)),
bottom= textGrob(expression(paste(log[10],'(',italic(mass),')',' [g]')), gp = gpar(cex = 1.8)))
# Combine and annotate a collection of the plot elements

dev.new(width=18,height=9,unit=cm)
ggarrange(h, annotated, labels=c('a','b'),font.label=list(size=25))
# Produce a complete combined plot

# Set the work directory to the location where you wish to save the combined plots.
ggsave(filename='diagnostic_plot_01_11_2023.pdf')
# Save the plots 
