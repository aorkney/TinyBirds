# This script runs analyses to generate Figure 2 and 3 from Orkney & Hedrick 2024.
# This code was written by Dr. Andrew Orkney, parsed by Prof. Brandon Hedrick, 
# and draws upon existing functions from dependent R packages and online contributions from Stackexchange.

# The purpose of this script is to load 
# centroid size information for 13 skeletal elements across a sample
# of 228 birds from across crown Aves. 
# The landmark constellations from which these data were extracted 
# are available at: https://doi.org/10.1038/s41586-022-05372-y
# (Navalon et al., 2022)
# Thereafter, the birds' masses will be used to perform
# a phylogenetically-informed Generalised Least Squares regression, 
# and the residuals will be extracted.
# This will remove allometric scaling in centroid size. 
# Residual centroid sizes of skeletal elements across birds are known 
# to exhibit a modular organisation of evolutionary covariances;
# some bones are more likely to evolve in concert than others. 
# This has widespread ramifications for how birds adapt to ecological demands. 
# https://doi.org/10.1038/s41559-021-01509-w
# (Orkney, Bjarnason, et al., 2021)
# This script will explore whether this modular organisation varies within birds, 
# specifically whether body mass variety is associated with different modular
# schemes of evolution. 

# Load required analytical R packages
library( geomorph ) # version 4.0.5
library( ape ) # version 5.0
library( nlme ) # version 3.1-162
library( phytools ) # version 1.5-1
# Done

setwd('C:/Users/Lab/Documents/Navalon_birds/File S1.Bird_landmarks_Navalon_2022/Clean')
# Set work directory appropriately (you will need to adjust this if you
# are replicating the work undertaken here)

# Load centroid sizes for study birds. 
load('Csize.12.4.2023.RData')
# In our study, we defined this object as an array called 'GPA.Csize'

# Load phylogenetic tree of birds
load('tree.12.4.2023.RData')
# Pruned tree produced as described in Navalon et al., 2022
# It is a merger of Prum and Oliveros
# The object is a phylogeny called 'pruned.tree'.

# Load bird masses
load('masses.12.4.2023.RData')
# The object is a named numeric vector called 'masses'.

# Define a sundry function to interract with data objects such as GPA.coords
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute the residual skeletal element centroid cizes, after the affects of allometry are removed
get.residual.Csize <- function( array, masses, phylogeny, taxa ){ # The function takes a shape array, mass vector, phylogeny and list of taxa
	species<-taxa
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to desired taxa
	for(i in 1:length(array) ){ # For each skeletal element 
		lambda <- phylosig(tree= newphy, x=log10( array[[ i ]][taxa]),  method='lambda')[[1]] # Initial estimate of phylogenetic autocorrelation
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  ) # Define a data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel(lambda, phy=newphy, form= ~ species ), data=df ) # Compute an allometric model
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

# This script will produce '.csv' files for 20 different bins of body mass, 
# showing how integration among pairwise combinations of bones changes with mass
# For the studied birds (228 species), these bins will be overlapping, 
# which will generate statistical non-independence between some bins
# This can be regarded as being similar to a rolling-average
# Because the birds are not uniformly distributed in frequency over log10(bodymass), 
# some bins will have wider or narrower ranges of mass represented than others

# Load package 'abind'; the 'array bind' function will be needed;
# If the user has not installed this package run 'install.packages('abind')' and proceed as guided
# no further comments advising package installation will be provided
library(abind) # version 1.4-5
# Package loaded

# Retrieve the names of the skeletal elements to be considered
bones<-names(residual.Csize)
# Names defined

# The following lines specify a function 
# that computes pairwise Z-scores (effect sizes) of integration between 
# different combinations of bones, for a stated group of taxa and a provided phylogeny relating the taxa

pair.int.phy <- function( input, phylogeny, taxa ){ # Define function names and required arguments; landmark coordinates (input), phylogeny and taxa

	newphy<-drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa] ) # Prune the source phylogeny as appropriate
	
	bones <- names(input) # Define the bone names in the input data
	output <- matrix(NA,length(bones),length(bones)) # Define a matrix to receive pairwise combinations of Z-scores between different bones
	output.p <- matrix(NA,length(bones),length(bones)) # Define  matrix to receive corresponding P-values
	rownames(output)<-bones # Name the rows and columns appropriately after the bones they represent
	colnames(output)<-rownames(output) # Set variable names (bone combinations) 
	colnames(output.p)<-colnames(output) # Set names for p-value matrix
	rownames(output.p)<-rownames(output) # Set names for p-value matrix
	bone_combinations<- combn(bones,2,simplify=T) # Define a vector of all possible combinations of bones that need to be considered
	
	int_list<-list() # Define a dumby list to receive integration estimates 
	for(i in 1:dim(bone_combinations)[2] ){ # For each pairwise combination 
		x<-paste(bone_combinations[,i][1]) Bone 1
		y<-paste(bone_combinations[,i][2]) Bone 2
		int_list[[i]] <- phylo.integration(as.matrix(input[[which(bones==x)]][newphy$tip]),as.matrix(input[[which(bones==y)]][newphy$tip]),phy=newphy,iter=999,print.progress=F) # Compute the Z-score and store it
	}

	#comparison<-compare.pls(int_list) # Compute the difference of integration statistics between the pairwise Z-scores
	
	z.score <- unlist( lapply(int_list,function(x) x$Z ) ) # Turn the resultant effect sizes of difference of integration into a vector

	p.val<-unlist( lapply(int_list,function(x) x$P.value ) ) # Extract the corresponding p-values 

	for(i in 1:dim(bone_combinations)[2] ){ # For each indiviual pairwise combination of bones
		x<-paste(bone_combinations[,i][1]) # Bone 1
		y<-paste(bone_combinations[,i][2]) # Bone 2
		output[which(rownames(output)==x),which(colnames(output)==y)]<-z.score[i] # Assign the relevant Z-score to the output matrix
		output[which(rownames(output)==y),which(colnames(output)==x)]<-z.score[i] # "
	}

	for(i in 1:dim(bone_combinations)[2] ){ # For each indiviual pairwise combination of bones
		x<-paste(bone_combinations[,i][1]) # Bone 1
		y<-paste(bone_combinations[,i][2]) # Bone 2
		output.p[which(rownames(output.p)==x),which(colnames(output.p)==y)]<-p.val[i] # Assign the corresponding P-value to the output matrix
		output.p[which(rownames(output.p)==y),which(colnames(output.p)==x)]<-p.val[i] # "
	}
		
	return(abind( output, output.p,along=3 )) # Return the Z-scores and corresponding p-values as a bound array
} # Function concludes

# Define the upper-bound of the body mass bins to be studied
breaks<-round(seq(40,length(masses),length.out=20))
# Bounds defined

# Define the width of the body mass bins to be studied, by number of taxa
bin.width<-40
# Width defined

# Define the number of taxa to be sampled in any individual replicate within a bin
n<-30 
subsample.size<-30
# Each replicate will be computed by subsampling 30 taxa within each bin

# Extract the number of pairwise combinations of skeletal elements
k <- dim(combn(names(residual.Csize),2,simplify=2))[2]
# Done

# Define matrices to receive the Z-scores, subsessetted by binned body mass, and associated p-values
size.scale.integration<-matrix( NA, length(breaks)*n, k+5 ) 
size.scale.integration.p<-matrix( NA, length(breaks)*n, k+5 ) 
# Matrices defined

# Define a vector that corresponds to the breaks between the body mass bins
break.vector<-rep(breaks,each=n)
# Done

# Comments will follow each line directly within the next loop, unless the function of plural lines is described, in which case the comment will be underneath the appropriate lines

for( i in 1 : dim(size.scale.integration)[1] ){ # For each replicate, within each body mass bin...

	end<-break.vector[i]
	start<-end-bin.width
	# Define the bin of body mass.

	bin.mass <- mean(masses[order(masses)][start:end])
	bin.birds <- names(masses[order(masses)][start:end])
	# Compute the mean mass and name the birds within the bin 

	selection<-sample(bin.birds,subsample.size)
	# Extract a subsample of taxa.

	newphy<-drop.tip( pruned.tree, pruned.tree$tip.label[ !pruned.tree$tip.label  %in% selection] ) # Prune the source phylogeny as appropriate
	phy.dist <- mean(cophenetic.phylo(newphy)) # Compute mean phylogenetic distance (we don't actually need this, but some users might appreciate it for diagnostics)
	
	subsample.mass <- mean(masses[order(masses)][start:end][selection])
	min.mass <- min(masses[order(masses)][start:end][selection])
	max.mass <- max(masses[order(masses)][start:end][selection])
	# Compute the mean, min and max body masses within the subsample of taxa. 

	size.scale.integration[i,1]<-bin.mass
	size.scale.integration[i,2]<-subsample.mass
	size.scale.integration.p[i,1]<-bin.mass
	size.scale.integration.p[i,2]<-subsample.mass
	size.scale.integration[i,3]<-min.mass
	size.scale.integration[i,4]<-max.mass
	size.scale.integration[i,5] <- phy.dist

	size.scale.integration.p[i,3]<-min.mass
	size.scale.integration.p[i,4]<-max.mass
	size.scale.integration.p[i,5] <- phy.dist
	# Store these information to the output matrix 

	result<-pair.int.phy( residual.Csize, pruned.tree, selection)# Compute pairwise integration between all combinations of skeletal elements

	size.scale.integration[i,-c(1:5)] <- result[,,1][upper.tri(result[,,1])]
	size.scale.integration.p[i,-c(1:5)] <- result[,,2][upper.tri(result[,,2])]
	# Store the resultant test statistics in the data output matrix

	print(round(i/dim(size.scale.integration)[1]*100,digit=2))# Print the progress to the screen (this code will take some time to run)
}
# This loop has produced a matrix of Z-scores and p-values that change with increasing body mass of the taxanomic subsets

# Define the index for the upper triangle of pairwise integration scores
ind<- which( upper.tri(result[,,1],diag=F) , arr.ind = TRUE )
# The lower triangle is identical, so is not needed 

# Define column names for the data output matrix
colpaste<-paste(bones[ind[,1]],bones[ind[,2]],sep='.')
# Done

# Set the column names for the objects recording Z-scores and p-values
colnames(size.scale.integration) <- c('bin.mass','subsample.mass','min.mass','max.mass','phy.dist',colpaste)
colnames(size.scale.integration.p) <- c('bin.mass','subsample.mass','min.mass','max.mass','phy.dist',colpaste)
# Assigned

# Check the working directory is as expected
getwd()
# User action may be warranted; this is where the files are to be saved 
setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_1_09_2024')

# Save the output for subsequent analyses
write.csv( x=size.scale.integration, file='size_scale_integration_01_15_2024_fixed.csv')
write.csv( x=size.scale.integration.p, file='size_scale_integration_p_01_15_2024_fixed.csv')
# argument 'x' describes the object, 'file' is the resultant file name

# Load packages that are needed to organise and plot data;
# If the user has not already installed a required package then run 
# 'install.packages('package name') and proceed as guided
library(ggplot2) # version 3.4.1 essential plotting functions
library(reshape2) # version 1.4.4 data organisation
library(ggpubr) # version 0.6.0 functions to combine plots
library(ggnewscale) # version 0.4.9 functions to overlay multiple fields of continuous data in plots
# Packages loaded

# Read the .csv file describing skeletal integration and its evolution with body mass
size.int<-read.csv('size_scale_integration_01_15_2024_fixed.csv')
# File read

# Remove the first row and set negative effect sizes to zero
size.int<-size.int[,-1]
size.int[size.int < 0 ]<-0
size.int[,-c(1:4)]<-size.int[,-c(1:4)]/(30^.5)
# Done

# Read the .csv file describing the P-values associated with skeletal integration and its evolution with body mass
size.int.p<-read.csv('size_scale_integration_p_01_15_2024_fixed.csv')
# File read

# Remove the first row and setvalues to zero
size.int.p<-size.int.p[,-1]
size.int.p[size.int.p < 0 ]<-0
# Done 

# Optional: let's try setting 'non-significant' values to zero
# size.int[,-c(1:4)][size.int.p[,-c(1:4)] >0.05] <- 0 

# The following lines define a vector of bones belonging to the 'head' module, 
# and retrieve the required indices to access the integration information
head.bones<-
combn(c('skull','mandible'),2) 
within.head<-matrix(NA,1,dim(head.bones)[2])
for(i in 1:dim(head.bones)[2] ){
	within.head[,i]<-intersect(grep(head.bones[1,i],colnames(size.int)),
	grep(head.bones[2,i],colnames(size.int)) )
}
within.head<-as.vector(within.head)
# Done

# The following lines define a vector of bones belonging to the 'wing' module, 
# and retrieve the required indices to access the integration information
wing.bones<-
combn(c('humerus','radius','ulna','carpometacarpus'),2) 
within.wing<-matrix(NA,1,dim(wing.bones)[2])
for(i in 1:dim(wing.bones)[2] ){
	within.wing[,i]<-intersect(grep(wing.bones[1,i],colnames(size.int)),
	grep(wing.bones[2,i],colnames(size.int)) )
}
within.wing<-as.vector(within.wing)
# Done

# The following lines define a vector of bones belonging to the 'trunk' module (body), 
# and retrieve the required indices to access the integration information
trunk.bones<-
combn(c('scapula','coracoid','sternum','synsacrum'),2) 
within.trunk<-matrix(NA,1,dim(trunk.bones)[2])
for(i in 1:dim(trunk.bones)[2] ){
	within.trunk[,i]<-intersect(grep(trunk.bones[1,i],colnames(size.int)),
	grep(trunk.bones[2,i],colnames(size.int)) )
}
within.trunk<-as.vector(within.trunk)
# Done

# The following lines define a vector of bones belonging to the 'leg' module, 
# and retrieve the required indices to access the integration information
leg.bones<-
combn(c('femur','tibiotarsus','tarsometatarsus'),2) 
within.leg<-matrix(NA,1,dim(leg.bones)[2])
for(i in 1:dim(leg.bones)[2] ){
	within.leg[,i]<-intersect(grep(leg.bones[1,i],colnames(size.int)),
	grep(leg.bones[2,i],colnames(size.int)) )
}
within.leg<-as.vector(within.leg)
# Done


# The following lines define a matrix to receive integration values binned by body
# module and body mass
inc<-unique(size.int$bin.mass)
mean.size.int<-matrix(NA,length(inc),21)
mean.size.int[,1]<-log10(inc)
colnames(mean.size.int)<-c('bin.mass',
'head','headmin','headmax','headminmin','headmaxmax',
'wing','wingmin','wingmax','wingminmin','wingmaxmax',
'trunk','trunkmin','trunkmax','trunkminmin','trunkmaxmax',
'leg','legmin','legmax','legminmin','legmaxmax')
# Done

# The commenting style will be simplified within the following loop, 
# which will produce boot-strapped probability intervals for change in integration within
# modular body regions, over body mass

for(i in 1: length(inc) ){ # For all unique mean body mass values
	ind<-which(size.int$bin.mass==inc[i]) # Find the associated indices
	mean.size.int[i,2]<-mean( size.int[ind,within.head] ) # Compute mean within-head integration

	boot<-list() # Define a dumby variable
	for(j in 1:1000){ # For 1000 replicates
		 boot[[j]]<-sample(size.int[ind,within.head],replace=T) # Sample integration with replacement
	}
	boot<-unlist(boot)
	mean.size.int[i,3]<-quantile(boot,probs=0.32) # Use the quantiles to extract probability intervals at 1 and 2 standard deviations
	mean.size.int[i,4]<-quantile(boot,probs=0.68)
	mean.size.int[i,5]<-quantile(boot,probs=0.05)
	mean.size.int[i,6]<-quantile(boot,probs=0.95)
	# Comments will not be repeated unnecessarily

	# Repeat process for within-wing integration:
	mean.size.int[i,7]<-mean( rowMeans( size.int[ind,within.wing] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.wing]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,8]<-quantile(boot,probs=0.32)
	mean.size.int[i,9]<-quantile(boot,probs=0.68)
	mean.size.int[i,10]<-quantile(boot,probs=0.05)
	mean.size.int[i,11]<-quantile(boot,probs=0.95)

	# Repeat process for within-trunk integration:
	mean.size.int[i,12]<-mean( rowMeans( size.int[ind,within.trunk] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.trunk]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,13]<-quantile(boot,probs=0.32)
	mean.size.int[i,14]<-quantile(boot,probs=0.68)
	mean.size.int[i,15]<-quantile(boot,probs=0.05)
	mean.size.int[i,16]<-quantile(boot,probs=0.95)

	# Repeat process for within-leg integration:
	mean.size.int[i,17]<-mean( rowMeans( size.int[ind,within.leg] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample( rowMeans(size.int[ind,within.leg]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,18]<-quantile(boot,probs=0.32)
	mean.size.int[i,19]<-quantile(boot,probs=0.68)
	mean.size.int[i,20]<-quantile(boot,probs=0.05)
	mean.size.int[i,21]<-quantile(boot,probs=0.95)

} # Loop concludes
# The output matrix has now been populated with mean within-module integration statistics, 
# as well as 2 and 1 sigma confidence envelope statistics


# The following lines will re-organise the data for plotting
df<-as.data.frame(mean.size.int)
molten.df<-melt(df,id='bin.mass')
# Done

# Define a dataframe for the confidence envelopes
ribbons<-c(grep('min',molten.df$variable),grep('max',molten.df$variable))
# Done

# Incorporate confidence envelopes into the re-organised data
melt.df<-melt(df[-ribbons,],id='bin.mass')
# Done 

# Define a colour blind friendly palette for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Done

mp<- # This plot will be produced so that its legend can be extracted. It is effectively a diagnostic plot
ggplot()+ # Make plot 
geom_line(aes(x=bin.mass,y=value,group=variable,linetype=variable,col=variable),size=1,data=molten.df[-ribbons,])+ # Plot lines of average within-module integration over body mass
scale_colour_manual(values=c('head'=cbbPalette[7],'wing'=cbbPalette[4],
'trunk'=cbbPalette[6],'leg'=cbbPalette[8]))+ # Define line colours of body modules 
labs(colour='module',linetype='module')+ # Label the line colours by corresponding body modules
# Sundry thematic settings follow that regulat plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(3, "cm"),
	legend.key.height = unit(1.5/sqrt(2), "cm"),
	legend.text=element_text(size=20),
	legend.title=element_text(size=20))
# A line plot has been produced (the purpose of this plot is just to extract its legend)

# Extract the plot's legend
mp.legend<-cowplot::get_legend(mp)
# Make sure the cowplot package is installed. 

# The first subplot will be commented for illustrative purposes; after this further subplot code will be identified at its start, but not commented explicitly
wing.mp<- # Figure 1 subplot a; evolution of within-wing integration as a function of body mass
ggplot()+ # Make plot 
geom_ribbon(aes(x=bin.mass,ymin=wingminmin,ymax=wingmaxmax),alpha=1,data=df,bg="#80FEF3")+ # Plot the 2 sigma confidence envelope of within-wing integration with body mass
geom_ribbon(aes(x=bin.mass,ymin=wingmin,ymax=wingmax),alpha=1,data=df,bg="#40CEB3")+ # Plot the 1 sigma confidence envelope of within-wing integration with body mass
geom_line(aes(x=bin.mass,y=wing),size=1,col=cbbPalette[4], data=df,linetype='22')+ # Plot a line tracking the mean within-wing integration with body mass
geom_text(x=2.0, y=5.25, aes(label=c('Wing')), size=6 )+ # Add a text label to represent the body module being considered; the wing
#lims(y=c(2.5,5.5))+ # Define axis limits
lims(y=c(0.35,0.95))+
labs(colour='module',linetype='module')+ # This line is not required; it was included as a result of re-using code, but it does not have any undesireable effect
labs(x='',y='')+
scale_x_continuous(position = "top")+
#labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+ # Label the x and y axes of the plot (currently disabled)
# Sundry instructions follow for the plot's thematic appearance:
  theme(axis.line.x.top=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plot produced

df$headmaxmax[which(df$headmaxmax>1)]<-1

head.mp<- # Figure 1 subplot d will be produced; it represents the evolution of within-head integration as a function of body mass
ggplot()+
#coord_fixed(ratio=(2/3)/(1.75/3))+ # Custom aspect ratios may be set if the user desires
geom_ribbon(aes(x=bin.mass,ymin=headminmin,ymax=headmaxmax),alpha=1,data=df,bg="#F5CE80")+
geom_ribbon(aes(x=bin.mass,ymin=headmin,ymax=headmax),alpha=1,data=df,bg="#F58E40")+
geom_line(aes(x=bin.mass,y=head),size=1,col=cbbPalette[7], data=df,linetype='solid')+
#lims(y=c(4.25,6.0))+ # Custom axis limits may be set if the user desires
geom_text(x=2.0, y=5.7, aes(label=c('Head')), size=6 )+
labs(colour='module',linetype='module')+
labs(x='',y='')+
scale_y_continuous(position = "right", limits=c(0.35,1))+
scale_x_continuous(position = "top")+
#labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+ # axis labels currently disabled 
  theme(axis.line.x.top=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.top=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

trunk.mp<- # Figure 1 subplot b will be produced; it represents the evolution of within-trunk integration as a function of body mass
ggplot()+
#coord_fixed(ratio=(2/3)/(3.25/3))+
geom_ribbon(aes(x=bin.mass,ymin=trunkminmin,ymax=trunkmaxmax),alpha=1,data=df,bg="#80E2F2")+
geom_ribbon(aes(x=bin.mass,ymin=trunkmin,ymax=trunkmax),alpha=1,data=df,bg="#40A2F2")+
geom_line(aes(x=bin.mass,y=trunk),size=1,col=cbbPalette[6], data=df,linetype='31')+
#lims(y=c(0.75,4))+
#lims(y=c(0.15,1))+
lims(y=c(0.2,1))+
geom_text(x=2.0, y=3.45, aes(label=c('Trunk')), size=6 )+
labs(colour='module',linetype='module')+
labs(x='',y='')+
scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits=c(0.16,1))+
#labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())



leg.mp<- # Figure 1 subplot e will be produced; it represents the evolution of within-leg integration as a function of body mass
ggplot()+
#coord_fixed(ratio=(2/3)/(2.25/3))+
geom_ribbon(aes(x=bin.mass,ymin=legminmin,ymax=legmaxmax),alpha=1,data=df,bg="#FCE9F7")+
geom_ribbon(aes(x=bin.mass,ymin=legmin,ymax=legmax),alpha=1,data=df,bg="#FCA9E7")+
geom_line(aes(x=bin.mass,y=leg),size=1,col=cbbPalette[8], data=df,linetype='33')+
#lims(y=c(2.5,4.75))+
scale_y_continuous(position = "right", limits=c(0.35,1))+
geom_text(x=2.0, y=4.4, aes(label=c('Leg')), size=6 )+
labs(colour='module',linetype='module')+
labs(x='',y='')+
#labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


# The following lines will assemble all the subplots into one combined plot 
# Define a left column of subplots
mp.left<-ggarrange(wing.mp,trunk.mp,ncol=1,nrow=2,
labels=c('a','b'),label.x=0,label.y=1,font.label=list(size=25))
blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done

library(grid)# base


# Define a right column of subplots
mp.right<-ggarrange(head.mp,leg.mp,ncol=1,nrow=2,
labels=c('c','d'),label.x=0.945,label.y=1,font.label=list(size=25))
# Done


# We are going to super-impose these plots upon 
# an illustrated background image

library(jpeg)# version 0.1-10
 
img <- readJPEG("bird_pastey.jpg")
# The user will not possess this image; they can dispense with running these lines

library(cowplot) # version 1.1.1

blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done

library(magick) # version 2.8.2

plot.with.inset <-
  ggdraw() +
  draw_plot(blank) +
  draw_image(img,x=0.22,y=0.21,width=0.57,height=0.57)+ # The user seeking to replicate our results should redact this line; it simply adds an illustration to the background
  draw_plot(mp.left, x = 0.01, y = 0, width = .48, height = 1)+
  draw_plot(mp.right, x = 0.50, y = 0, width = .48, height = 1.005)+
	draw_plot(mp.legend,x=0.11,y=0.38,width=0.2,height=0.3)


annotation.with.inset <- annotate_figure(plot.with.inset, left = textGrob( expression(paste( italic('Z/'),sqrt(n) ) ), rot = 90, vjust = 0.5,  gp = gpar(fontsize=20)),
bottom= textGrob(expression(paste(log[10],'(',italic(mass),')',' [g]')), gp = gpar(fontsize=20)))
# Label the plot

dev.new(height=10,width=12.25,unit='cm')
annotation.with.inset
# Produce the final plot

# Save the plot as a pdf and as a high quality .PNG file
ggsave(filename='mass_v_integration_01_11_2024.pdf')
# Done




# The following lines of code will compute Ordinary Least Squares fits for change in pairwise-Z scores
# of integration across all combinations of skeletal elements as a function of mean body mass
# within each bin
# This will produce Figure 3s, which is a plot produced for illustrative purposes; be aware that the bins of body mass
# may overlap, causing statistical non-independence among the data

k<-78
# 78 combinations of bones
kseq<-inc[seq(1,20,by=1)]
# A sequence of unique body mass values 

# These lines define a matrix to receive the slopes of OLS fits:
recompile <- matrix(NA,13,13)
rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
colnames(recompile)<-rownames(recompile)
# Done

# The following loop will populate the matrix with OLS slope fits:
inc<-unique(size.int[,1])
for(i in 1:length(size.int[1,-c(1:2)])  ){
	names<-colnames(size.int[,-c(1:2)])
	x<-  sub("\\..*","", names)[i]
	y<-  sub(".*\\.","", names)[i]
	ind<-which(size.int$bin.mass==kseq[k])

	temp<-list()
	for(j in 1:length(inc)){
		ind<-which(size.int[,1]==inc[j])
		temp[[j]]<-mean(size.int[,-c(1:2)][,i][ind])
		
	}
	temp<-unlist(temp)

	coef<-lm(temp~log10(inc))

	if(summary(coef)$coef[2,4] <1){ # The user may exclude cells above a specific significance threshold if they wish, by replacing '<1'
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-coef$coef[2]
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-coef$coef[2]
	} 
}
# Done

# Re-organise the output for plotting
result.rev<-apply(t(recompile),2,rev)
longData<-melt(result.rev)
# Done

levels(longData$Var1)[13]<-'cranium'
levels(longData$Var2)[1]<-'cranium'
longData$Var1[is.na(longData$Var1)]<-'cranium'
longData$Var2[is.na(longData$Var2)]<-'cranium'
# Changing the data column named 'Skull' to 'cranium'


# Define colours for different body modules of the avian skeleton 
bone_colours<-rev(c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3)))
# Done

mp.change<- # Produce Figure 1 subplot c
ggplot(longData, aes(x = Var2, y = Var1)) + # Make plot, specify data sources and axial variables
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ # Define geometry 
  geom_tile(aes(fill=value)) + # Treat the plot as a tile plot 
   scale_fill_gradient2(low='blue',mid='white',high='red' ,na.value="white",midpoint=0) + # Fill the tiles with colour corresponding to the slope of OLS fit
	labs(fill = "slope")+ # Label the legend
  labs(x="", y="") + # No axial labels are needed
# Sundry thematic appearance instructions follow
#scale_x_discrete(guide = guide_axis(angle = 90), size=14)+ # Custom x axis option
  theme_bw() + theme(axis.text.x=element_text(angle=90, size=14,hjust=0, colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(size=14,colour=(bone_colours),face = "bold"),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=20),
				legend.text=element_text(size=14),
				plot.margin=unit(c(-0.5,-1,0.5,-1),'cm'),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = expression(paste(Delta,italic('Z/'),sqrt(n) )))+
geom_hline(yintercept = 3.5,size=3/2,colour="black")+
geom_hline(yintercept = 11.5,size=3/2,colour="black")+
geom_hline(yintercept = 7.5,size=3/2,colour="black")+
geom_vline(xintercept = 10.5,size=3/2,colour="black")+
geom_line(data=as.data.frame(cbind(c(2.5,2.5),c(0.5,13.5))),aes(x=V1,y=V2),size=2,colour='black')+
geom_line(data=as.data.frame(cbind(c(6.5,6.5),c(0.5,13.5))),aes(x=V1,y=V2),size=2,colour='black')+
geom_line(data=as.data.frame(cbind(c(10.5,10.5),c(0.5,13.5))),aes(x=V1,y=V2),size=2,colour='black')+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=2.5,ymin=11.5,ymax=13.5,linewidth=3.5)+
annotate('rect',colour='#D55E00',fill=NA,xmin=0.5,xmax=2.5,ymin=11.5,ymax=13.5,linewidth=3/2)+
annotate('rect',colour='black',fill=NA,xmin=2.5,xmax=6.5,ymin=7.5,ymax=11.5,linewidth=3.5)+
annotate('rect',colour='#009E73',fill=NA,xmin=2.5,xmax=6.5,ymin=7.5,ymax=11.5,linewidth=3/2)+
annotate('rect',colour='black',fill=NA,xmin=6.5,xmax=10.5,ymin=3.5,ymax=7.5,linewidth=3.5)+
annotate('rect',colour='#56B4E9',fill=NA,xmin=6.5,xmax=10.5,ymin=3.5,ymax=7.5,linewidth=3/2)+
annotate('rect',colour='black',fill=NA,xmin=10.5,xmax=13.5,ymin=0.5,ymax=3.5,linewidth=3.5)+
annotate('rect',colour='#CC79A7',fill=NA,xmin=10.5,xmax=13.5,ymin=0.5,ymax=3.5,linewidth=3/2)+
scale_x_discrete(position = "top") 
# Plot complete
# This has produced a plot that shows general trends in pairwise-Z of integration across the skeleton as a function
# of increasing body mass

# It is time to make microplots which demonstrate how pairwise schemes of integration change across different body mass ranges, with a collection of 'snapshots'. 

# The following 'yvals' define the arbitrary axial position of the microplots
# Here we are going to line them all up along the same value
yvals<-rep(1,6)

# 'kseq' is the vector of microplot identities that will be plotted, corresponding to x-axial positions
# The user may make bespoke choices about which microplots to represent, but I have selected an aesthetic and roughly evenly spaced arrangement:
kseq<-inc[c(3,9,14,18)]


blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done

mp.comp<-blank+
#scale_y_continuous(position = "right")+
labs(y=expression(paste(log[10],'(',italic(mass),')',' [g]')), x='')+
scale_y_continuous(breaks=log10(kseq),
        labels=c('1.4','1.9','2.5','3.0'), position = "right")+
  theme(
	legend.key.width=unit(1/2,'cm'),
	legend.key.height=unit(1,'cm'), 
	legend.title=element_text(size=20),
	legend.text=element_text(size=14),
	axis.line.x.top=element_line(size = 1, colour = "black"),
	axis.line.y.right=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x=element_blank(),
	axis.title.y.right=element_text(size=20,colour="black", vjust = 2),
      axis.ticks.x=element_blank(),
	axis.ticks.y=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

for(k in 1:length(kseq)){ # For each microplot
	recompile <- matrix(NA,13,13) # The microplot will be 13 by 13 cells wide
	rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
	'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
	colnames(recompile)<-rownames(recompile) # The cells in the microplot correspond to average integration between pairwise combinations of bones across the bird body plan

	for(i in 1:length(size.int[1,-c(1:2)])  ){ # This loop will extract the relevant integration data from a previously produced and saved .csv file
		names<-colnames(size.int[,-c(1:2)])
		x<-  sub("\\..*","", names)[i]
		y<-  sub(".*\\.","", names)[i]
		ind<-which(size.int$bin.mass==kseq[k])
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-mean( size.int[ind,-c(1:2)][,i] )
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-mean( size.int[ind,-c(1:2)][,i] )
	} # Loop concludes

	result.rev<-apply(t(recompile),2,rev) # The resultant data for the microplot will be transposed
	longData<-melt(result.rev) # The data will be re-organised into a format suitable to supply to plotting functions
	longData$value[which(longData$value<0)]<-0 # Cells which have a value below 0 will be set to 0 
	longData$Var1<-(as.numeric(longData$Var1)/(100/(5.5-2.5)))+log10(kseq[k]) # The y-axial position of the microplot will be specified
	longData$Var2<-(as.numeric(longData$Var2)/(100/(2))) +yvals[k] # The x-axial position of the microplot will be specified

	scale.y <- diff(range(longData$Var1))
	scale.x <- diff(range(longData$Var2))

	mp.comp<- mp.comp+ # The overlaying operation will be performed
	coord_fixed(ratio=2/3)+ # A set coordinate ratio is defined to conserve geometry
  	geom_raster(aes(fill=value,x = Var2, y = Var1),data=longData)+ # The overlaid microplot is defined as a raster object
	scale_fill_gradientn(colours = c('white','black') ,na.value="white",lim=c(0,1))+ # Stronger integration within the microplot will be represented by a darker colour
	labs(fill = expression(paste(italic('Z/'),sqrt(n) )))+
	new_scale("fill")+  # The new plotting scale of the overlaid microplot is defined as a colour fill
	annotate('rect',colour='black',fill=NA,xmin=1.01,xmax=1.051,ymin= unique(longData$Var1)[11]+0.014, ymax= unique(longData$Var1)[13]+0.015,linewidth=1.5)+
	annotate('rect',colour='#D55E00',fill=NA,xmin=1.01,xmax=1.051,ymin= unique(longData$Var1)[11]+0.014, ymax= unique(longData$Var1)[13]+0.015,linewidth=1)+
	annotate('rect',colour='black',fill=NA,xmin=1.05,xmax=1.131,ymin= unique(longData$Var1)[7]+0.014, ymax= unique(longData$Var1)[11]+0.015,linewidth=1.5)+
	annotate('rect',colour='#009E73',fill=NA,xmin=1.05,xmax=1.131,ymin= unique(longData$Var1)[7]+0.014, ymax= unique(longData$Var1)[11]+0.015,linewidth=1)+
	annotate('rect',colour='black',fill=NA,xmin=1.13,xmax=1.211,ymin= unique(longData$Var1)[3]+0.014, ymax= unique(longData$Var1)[7]+0.015,linewidth=1.5)+
	annotate('rect',colour='#56B4E9',fill=NA,xmin=1.13,xmax=1.211,ymin= unique(longData$Var1)[3]+0.014, ymax= unique(longData$Var1)[7]+0.015,linewidth=1)+
	annotate('rect',colour='black',fill=NA,xmin=1.21,xmax=1.271,ymin= log10(kseq[k])+0.014, ymax= unique(longData$Var1)[3]+0.015,linewidth=1.5)+
	annotate('rect',colour='#CC79A7',fill=NA,xmin=1.21,xmax=1.271,ymin= log10(kseq[k])+0.014, ymax= unique(longData$Var1)[3]+0.015,linewidth=1)



}# Loop condludes
mp.comp

dev.new(height=8,width=12,unit='cm')
ggarrange(blank,mp.change,mp.comp,ncol=3,widths=c(0.05,1,.5),labels=c('','a','b'),label.x=0,label.y=1,font.label=list(size=25))
 ggsave(filename='microplots_01_15_2024.pdf')

# Figure 3 produced and saved

