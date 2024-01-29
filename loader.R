# This script was originally written by Dr Roger Benson and has been modified by Dr Andrew Orkney

setwd( "C:/Users/Lab/Documents/Navalon_birds/File S1.Bird_landmarks_Navalon_2022" )	##Set as working directory the "Bird_landmarks_Navalon_2022" folder (User may need to customise)
	
library(ape)
library(geomorph)

tree <- read.tree('Navalon_tree.tre')
load('Navalon_birds.RData')

SL.counts <- "min"

elements <- c( "skull" , "mandible" , "scapula", "coracoid", "sternum" , "humerus", "ulna", "radius", "carpometacarpus", "synsacrum" , "femur", "tibiotarsus", "tarsometatarsus" )

name_list<-list()
for(i in 1:length(elements)){
	name_list[[i]]<-dimnames(compiled.bird.landmarks.list[[ elements[ i ] ]][[ SL.counts ]][[ "output.landmarks.array" ]] )[[ 3 ]]
}


AA<-table(unlist(name_list))
birds <- names( AA )[ AA == length( elements ) ]

GPA.results <- list()
	for( i in 1:length( elements ) )	{
		element <- elements [ i ]
			GPA.results[[ i ]] <- gpagen( compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "output.landmarks.array" ]][ , , birds ] ,  curves = compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "sliders" ]] , ProcD = T )
		}
names( GPA.results ) <- elements


get.item <- function( X , item ) { X[[ item ]] }	
GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )

# We have organised the data into a single object

metadata <- read.csv('Navalon_metadata.csv')

masses<- rowMeans( metadata[,11:12])

names(masses) <- metadata$Data.name

masses<- masses[match(names(GPA.Csize$sternum),names(masses))]

# Now it is time to save the data
setwd('C:/Users/Lab/Documents/Navalon_birds/File S1.Bird_landmarks_Navalon_2022/Clean') # The user may need to customise

pruned.tree <- tree

save(masses, file='masses.12.4.2023.RData')
save(pruned.tree,file='tree.12.4.2023.RData')
save(GPA.Csize,file='Csize.12.4.2023.RData')
save(GPA.coords,file='coords.12.4.2023.RData')
