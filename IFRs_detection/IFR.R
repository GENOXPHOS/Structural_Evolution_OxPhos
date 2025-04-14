library("bio3d")
library("parallel")

# Define mtDNA-encoded and X-linked subunits.

MT.subunits=c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712", 
           "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", 
           "ENSG00000212907","ENSG00000198886", "ENSG00000198786", "ENSG00000198695", 
           "ENSG00000198727")

X.subunits=c("ENSG00000125356", "ENSG00000131174", "ENSG00000147123")


##################################################
### GETS IFR FROM IN SILICO MSTRUCTURES OF RCs ###
##################################################

# Define class for object S4 where all information will be stored
setClass("Interfaces", representation(name = "vector", Ens.subunit.query = "vector", Ens.subunit.target="vector", subunit.query = "vector",
                                      subunit.target="vector", Interface.residues.query="list",genome.query="vector", genome.target="vector"))

# Create object S4 of class interfaces
IFR=new("Interfaces",name=character(),Ens.subunit.query=character(),
       Ens.subunit.target=character(),subunit.query = character(),subunit.target=character(),
       Interface.residues.query=list(),genome.query=character(),genome.target=character())



# List structure files alongside their corresponding tsv files, which provide the mapping between each chain and its associated subunit
path="/path/to/input_directory"
l1=list.files(path,pattern = ".Chain.gene.key.tsv")
l2=list.files(path,pattern = ".pdb")

# Iterate files to be analyzed 
for (j in 1:length(l1)){
  
  # Load PDB structure along with its corresponding tsv file for mapping chains to subunit
  key=read.delim(paste0(path,"/",l1[j]),sep = "\t",stringsAsFactors = F)
  pdb = read.pdb(paste0(path,"/",l2[j]))
  
  # Get name of studied respiratory complex (RC)
  name=gsub(".pdb","",l2[j])
  
  # Define all possible combinations of pairs chains (query-target) in the structure to be evaluated in search for IFRs
  grid=expand.grid(unique(pdb$atom$chain),unique(pdb$atom$chain))
  grid=grid[which(grid$Var1!=grid$Var2),]
  
  # Annotate S4 object with RC information
  IFR@name=c(IFR@name,name)
  
  # Annotate the S4 object with the Ensembl ID and the symbol of the gene representing the evaluated chains
  IFR@Ens.subunit.query=c(IFR@Ens.subunit.query,sapply(grid[,1],function(x)return(key[key[["Chain"]]==x,"ENSEMBL_ID"])))
  IFR@Ens.subunit.target=c(IFR@Ens.subunit.target,sapply(grid[,2],function(x)return(key[key[["Chain"]]==x,"ENSEMBL_ID"])))
  IFR@subunit.query=c(IFR@subunit.query,sapply(grid[,1],function(x)return(key[key[["Chain"]]==x,"Gene"])))
  IFR@subunit.target=c(IFR@subunit.target,sapply(grid[,2],function(x)return(key[key[["Chain"]]==x,"Gene"])))
  
  # Annotate the IFRs between chains 
  IFR@Interface.residues.query=c(IFR@Interface.residues.query,mclapply(1:nrow(grid),function(i){binding.site(pdb, a.inds=atom.select(pdb, chain=as.character(grid[i,1])), b.inds=atom.select(pdb, chain=as.character(grid[i,2])), verbose = F)$resno},mc.cores = 11))
}

# Annotate IFRs between chains with their values of encoding chromosomes for involved chains at the binding site
IFR@genome.query=IFR@Ens.subunit.query
IFR@genome.query[IFR@genome.query %in% X.subunits]="X"
IFR@genome.query[IFR@genome.query %in% MT.subunits]="Mitochondrial"
IFR@genome.query[IFR@genome.query %in% c("X","Mitochondrial") == F]="Autosomal"

IFR@genome.target=IFR@Ens.subunit.target
IFR@genome.target[IFR@genome.target %in% X.subunits]="X"
IFR@genome.target[IFR@genome.target %in% MT.subunits]="Mitochondrial"
IFR@genome.target[IFR@genome.target %in% c("X","Mitochondrial") == F]="Autosomal"

# Eliminates all non-neighboring combinations between chains 
cond=unlist(lapply(IFR@Interface.residues.query,function(x){return(length(x))}))>0
IFR@Ens.subunit.query=IFR@Ens.subunit.query[cond]
IFR@Ens.subunit.target=IFR@Ens.subunit.target[cond]
IFR@subunit.query=IFR@subunit.query[cond]
IFR@subunit.target=IFR@subunit.target[cond]
IFR@genome.query=IFR@genome.query[cond]
IFR@genome.target=IFR@genome.target[cond]
IFR@Interface.residues.query=IFR@Interface.residues.query[cond]

# Save S4 object containing information for IFRs
saveRDS(IFR,"/path/to/output_directory/IFRs.rds")
