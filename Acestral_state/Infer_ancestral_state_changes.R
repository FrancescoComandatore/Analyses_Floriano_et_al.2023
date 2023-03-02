
# This script has been written by Francesco Comandatore (University of Milan)
# If you use this script please cite:

#"The origin and evolution of mitochondrial tropism in Midichloria bacteria"
#Anna Maria Floriano, Gherard Batisti Biffignandi, Michele Castelli, Emanuela Olivieri, Emanuela Clementi, Francesco Comandatore, Laura Rinaldi, Maxwell Opara, Olivier Plantard, Ana M. Palomar, Valérie Noël, Amrita Vijay, Nathan Lo, Benjamin L. Makepeace, Olivier Duron, Aaron Jex, Lionel Guy, Davide Sassera

##############
# How to RUN
##############

# Rscript Infer_ancestral_state_changes.R [table file] [tree] [probability parameter for the ancestral state recontruction, def: 0.7]

############
# The scope 
############

# the ancestral state of the feature is reconstructed on the phylogenetic tree and the nodes corresponding to variations from a state to another are identified

#############
# The INPUTS 
#############

# 1. a ROOTED phylogenetic tree in Newick format
# 2. a tab delimited table in which the first column report the a label for each feature (e.g. the gene name) and each other column report the state of that feature for the organims (see examples below)

# All the organisms present in the tree must be present in the table



# Example 1 (binary)

#Gene	Org1	Org2	Org3
#gene1	0	1	1
#gene2	1	1	0
#gene3	1	1	1

# Example 2 (SNPs)

#Gene	Org1	Org2	Org3
#gene1	A	T	T
#gene2	T	G	A
#gene3	T	T	A


##############
# The OUTPUTS 
##############

# two differents output files are produced:

# 1. the output table which includes all the information of the input table but with added a column reporting the node in which a feature changed (including the original and modified states)
# 2. the labelled phylogenetic tree: the phylogenetic tree with the nodes labelled as in the output table. Opening this file with a visulization software (e.g. Seaview) it will be easy possible to identify the nodes relative to the feature changes

library(ape) 
library(phytools)
require(methods)

# Get inputs names, first annotation tab and then the ROOTED tree

args <- commandArgs(trailingOnly = TRUE)

# Opens input files

if (args[3] == ""){args[3] = 0.7}

### OPEN INPUTS

ricodata <- read.table(args[1],check.names=F,quote="",header=T, sep="\t")
ricotree <- read.tree(args[2]) 

# Get node and tip names

ricotree2=makeNodeLabel(ricotree)
node_tip_names=c(ricotree2$tip.label, ricotree2$node.label)

# Write the annotated tree

write.tree(ricotree2, file=paste(args[2],".labellednodes.tree",sep=""))
treematrix<-ricotree$edge

## Work on the dataframe

#get tree organisms columns

row.names(ricodata)<-ricodata[,1]

ricodata_tree <- ricodata[,ricotree$tip.label]

# polish on pattern

# create a column for each pattern
String_pattern<-gsub(" ","",do.call("paste",ricodata_tree))
ricodata_tree_pat<-cbind.data.frame(ricodata_tree,String_pattern)

ricodata_all_pat<-cbind.data.frame(ricodata,String_pattern)

# remove duplicated rows
ricodata3=unique(ricodata_tree_pat)
row.names(ricodata3)<-ricodata3$String_pattern
ricodata3$String_pattern <- NULL

### Perform analyses

output=matrix(ncol=4,nrow=0)
colnames(output)<-c("Node","String_pattern","STATE_from", "STATE_to")

for (z in 1:length(row.names(as.matrix(ricodata3))))
{
	variability=length(unique(as.factor(as.matrix(ricodata3)[z,])))

	if (variability > 1)
	{

	res <- try(fitER <- rerootingMethod(ricotree, as.matrix(ricodata3)[z,], model = "ER"))

	if (!inherits(res, "try-error"))
	{	
		fitER2<-(fitER$marginal.anc>args[3])

		base=NULL
		base=list()
		base<-ricodata3[z,]
		colnames(base)<-1:length(ricodata3[z,])

		p=length(ricotree$tip.label)+1

			for (i in 1:length(fitER2[,1]))
			{

				    if (sum(fitER2[i,]) == 1)
    					{
						base[p]<-names(subset(fitER2[i,],fitER2[i,]==1))        
        					colnames(base)[p]<-p
        					p = p + 1
    					}

				    else
    					{
        					base[p]<-NA
        					colnames(base)[p]<-p
        					p = p + 1
    					}
			}
	
		Nodes_comp<-na.omit(as.data.frame(cbind(t(base[treematrix[,2]]),t(base[treematrix[,1]]))))

		Nodes_sel<-Nodes_comp[as.vector(Nodes_comp[,1]) != as.vector(Nodes_comp[,2]), ] 
		Nodes_Pat<-cbind(node_tip_names[as.numeric(row.names(Nodes_sel))],rep(colnames(Nodes_sel)[1],length(row.names(Nodes_sel))), as.matrix(Nodes_sel[,2]), as.matrix(Nodes_sel[,1]))

		colnames(Nodes_Pat)<-c("Node","String_pattern","STATE_from", "STATE_to")

		output= rbind(output, Nodes_Pat)

		}


		}

}


Merged_table <- merge(output, ricodata_all_pat, by = "String_pattern")
Merged_table$String_pattern <- NULL

write.table(Merged_table, file=paste(args[1],".anc_state_changes.tab",sep=""),sep="\t",row.names=F, quote=F)


