library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)
#install.packages("RColorBrewer")
library(ratematrix)
library(OUwie)
library(RRphylo)
#library(RColorBrewer)



get_best_fit_model_AICc <- function(tree, neo) {
  best_model<-list()
  # Ensure the names of the response data (neo) match the tip labels in the tree
  if (!all(names(neo) %in% tree$tip.label)) {
    stop("The names of the response vector (neo) must match the tip labels in the phylogenetic tree.")
  }
  
  # Define the models to test
  models_to_test <- c("BM", "lambda", "OU", "EB", "delta", "kappa", "rate_trend")
  
  # Initialize a list to store the fit results
  fit_results <- list()
  
  # Try fitting each model
  for (model in models_to_test) {
    # First attempt: fit the model normally
    fit <- fitContinuous(
      phy = tree,
      dat = neo,
      model = model
    )
    
    # Store the fit if successful
    if (!is.null(fit)) {
      fit_results[[model]] <- fit
    }
  }
  
  # Check if we have any successful fits
  if (length(fit_results) == 0) {
    stop("None of the models could be fitted successfully.")
  }
  
  # Extract log-likelihood values for each successfully fitted model
  logLik_values <- sapply(fit_results, function(fit) fit$opt$aicc, USE.NAMES = TRUE)
  
  # Get the model with the highest log-likelihood
  best_model$allmodels <- logLik_values
  best_model$bestfit <- names(which.min(logLik_values))
  print(best_model)
  
  return(best_model$bestfit)
}



Data <- read.csv("min10ruminants.csv")

tree <- read.tree("min1.nwk")
#for malignancy

Data <- subset(Data, Family == "Bovidae" | Family == "Cervidae" | Family == "Giraffidae" | Family == "Camelidae" | Family == "Ovidae" | Family == "Antilocapridae")


cutDataAll <- Data[, c(9,17,7 ,6), drop=FALSE]
cutDataAll$Species <- gsub(" ", "_", cutDataAll$Species) 




# for neo
#cutData <- filterData[, c(9,13), drop=FALSE]

includedSpecies <- cutDataAll$Species
pruned.tree <- drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree, pruned.tree$tip.label)
cutDataAll$Keep <- cutDataAll$Species %in% pruned.tree$tip.label
cutDataAll <- cutDataAll[!(cutDataAll$Keep == FALSE), ]
rownames(cutDataAll) <- cutDataAll$Species
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
x <- setNames(as.numeric(cutDataAll[, 2]), cutDataAll[, 1])


best_model<-get_best_fit_model_AICc(pruned.tree, x)

tree<-multi2di(tree,random=TRUE)
is.binary(tree)

name.check(pruned.tree, x)

#edgelabels()
N<-length(pruned.tree$tip.label)




fit<-fitContinuous(pruned.tree, x, model = "OU")

#pagelFit<- fitContinuous(pruned.tree, x, model="lambda")


# 
# pruned.tree <- rescale(pruned.tree, model = "lambda", lambda = pagelFit$opt$lambda)
# pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10

name.check(pruned.tree, x)



#pagel<-anc.ML(pruned.tree, x, model = "BM")

ou<-fastAnc(pruned.tree, x, model = "OU")

#pagel<-pagel$ace

# 
#  pagel<-vector(); 
#  N<-length(pruned.tree$tip); M<-pruned.tree$Nnode
#  
#   
#  # compute the PIC "root" state for each internal node
#    for(i in 1:M+N){
#      pagel[i-N]<-ace(x,multi2di(root(pruned.tree,node=i)),
#                  method="pic")$ace[1]
#     names(pagel)[i-N]<-i
#   }


obj <- contMap(pruned.tree, x, method = "user", anc.states = ou, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "Malignancy Prev OU")




x_binary <- as.integer(x > 0)
names(x_binary)<-names(x)


pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10

#threshOU<-ancThresh(pruned.tree,x_binary, model = "lambda")

cutDataAll$Species <- gsub(" ", "_", cutDataAll$Species) 


# Filter your dataset by each family
bovidae <- filter(cutDataAll, Family == "Bovidae")
cervidae <- filter(cutDataAll, Family == "Cervidae")
giraffidae <- filter(cutDataAll, Family == "Giraffidae")
camelidae <- filter(cutDataAll, Family == "Camelidae")
ovidae <- filter(cutDataAll, Family == "Ovidae")
antilocapridae <- filter(cutDataAll, Family == "Antilocapridae")

# Set neutral colors for the labels
neutral_colors <- c(
  "Bovidae" = "#7F7F7F",          # medium gray
  "Cervidae" = "#A9A9A9",         # dark gray
  "Giraffidae" = "#BEBEBE",       # light gray
  "Camelidae" = "#D8BFAA",        # beige/taupe
  "Ovidae" = "#C0C0C0",           # silver gray
  "Antilocapridae" = "#E0D5C6"    # light taupe
)

tree_species <- pruned.tree$tip.label

bovidae_in_tree <- bovidae$Species[bovidae$Species %in% tree_species]


# Plot clade labels by family with neutral colors
par(fg=neutral_colors["Bovidae"])
cladelabels(pruned.tree,
  text = "Bovidae",
  cex = 1.05,
  node = findMRCA(pruned.tree, as.vector(bovidae$Species))
)

par(fg=neutral_colors["Cervidae"])
cladelabels(
  text = "Cervidae",
  cex = 1.05,
  node = findMRCA(pruned.tree, cervidae$Species)
)

par(fg=neutral_colors["Giraffidae"])
cladelabels(
  text = "Giraffidae",
  cex = 1.05,
  node = findMRCA(pruned.tree, giraffidae$Species)
)

par(fg=neutral_colors["Camelidae"])
cladelabels(
  text = "Camelidae",
  cex = 1.05,
  node = findMRCA(pruned.tree, camelidae$Species)
)

par(fg=neutral_colors["Ovidae"])
cladelabels(
  text = "Ovidae",
  cex = 1.05,
  node = findMRCA(pruned.tree, ovidae$Species)
)

par(fg=neutral_colors["Antilocapridae"])
cladelabels(
  text = "Antilocapridae",
  cex = 1.05,
  node = findMRCA(pruned.tree, antilocapridae$Species)
)

# Optionally, reset fg color to default after plotting
par(fg="black")
