rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./pop-dynamics-forams/analysis")

# Libraries
# source("R/library.R")
library(ape)
library(ggtree)
library(paleoPhylo) 

# Auxiliary functions
sourceDirectory("./R/aux_functions", modifiedOnly=FALSE)

### AZE

load(file = "data/phylogeny/2011-04-11aM.Rdata", envir = parent.frame(), verbose = FALSE)
names(aM) # name parent start end

load(file = "data/phylogeny/2011-04-11aL.Rdata", envir = parent.frame(), verbose = FALSE)
names(aL)

aze_tree <- buildApe(createBifurcate(aM))
aze_tree$tip.label[order(aze_tree$tip.label)]


write.tree(aze_tree, file = "data/phylogeny/Aze_2011-04-11aM.tre")


### KUCERA

# reading the text file with the raw tree (without node dates / branching times)
tree <- read.tree(file="data/phylogeny/phylogeny_21May.txt")

# plotting tree
ggtree(tree, layout ='rectangular', ladderize=TRUE) + geom_tiplab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, color = "blue") +
  xlim(0,15)

# saving the tree.txt in other formats (.nex or .tre)
# write.nexus(tree, file = "data/phylogeny/kucera_tree.nex")
# write.tree(plankforams, file = "data/phylogeny/kucera_tree.tre")

tree <- compute.brlen(tree, 1) # sets all branch lengths equal to one

ssp_dist <- cophenetic(tree) # returns a symmetrical matrix with distances for each species pair


dist.nodes(tree)
branching.times(tree)


polyroot <- read.tree(text="((spinose,nonspinose),(H._pelagica,H._digitata),micro,B._pumilio,D._anfracta);")
micro <- read.tree(text = "((((C._nitida,G._glutinata),G._uvula),G._minuta),(T._fleisheri,T._iota,T._parkerae));")
nonspinose <- read.tree(text ="(((((G._tumida,G._ungulata),G._menardii),((G._cavernula,G._truncatulinoides),G._crassaformis),(G._scitula,(G._theyeri,G._hirsuta))),G._inflata,((P._obliquiloculata,N._dutertrei),(N._incompta,N._pachyderma))),G._hexagonus, G._conglomerata);")
spinose <- read.tree(text ="((T._clarkei,(T._quinqueloba,T._humilis)),(((G._bulloides,G._falconensis),(B._digitata,(G._radians,G._calida,G._siphonifera,G._adamsi))),(((G._sacculifer,O._universa),S._dehiscens),(G._rubescens,(((G._elongatus,G._tenellus),G._conglobatus),(G._white,G._ruber))))));")
plot(spinose)


plankforams_kucera_original <- read.tree(text="((((Turborotalita_clarkei,(Turborotalita_quinqueloba,Turborotalita_humilis)),(((Globigerina_bulloides,Globigerina_falconensis),(Beella_digitata,(Globigerinella_radians,Globigerinella_calida,Globigerinella_siphonifera,Globigerinella_adamsi))),(((Trilobatus_sacculifer,Orbulina_universa),Sphaeroidinella_dehiscens),(Globoturborotalita_rubescens,(((Globigerinoides_elongatus,Globigerinoides_tenellus),Globigerinoides_conglobatus),(Globigerinoides_white,Globigerinoides_ruber)))))),(((((Globorotalia_tumida,Globorotalia_ungulata),Globorotalia_menardii),((Globorotalia_cavernula,Globorotalia_truncatulinoides),Globorotalia_crassaformis),(Globorotalia_scitula,(Globorotalia_theyeri,Globorotalia_hirsuta))),Globoconella_inflata,((Pulleniatina_obliquiloculata,Neogloboquadrina_dutertrei),(Neogloboquadrina_incompta,Neogloboquadrina_pachyderma))),Globorotaloides_hexagonus, Globoquadrina_conglomerata)),(Hastigerina_pelagica,Hastigerinella_digitata),((((Candeina_nitida,Globigerinita_glutinata),Globigerinita_uvula),Globigerinita_minuta),(Tenuitella_fleisheri,Tenuitella_iota,Tenuitella_parkerae)),Berggrenia_pumilio,Dentigloborotalia_anfracta);")

plankforams <- read.tree(text="(((((Turborotalita_clarkei,(Turborotalita_quinqueloba,Turborotalita_humilis)),(((Globigerina_bulloides,Globigerina_falconensis),(Beella_digitata,(Globigerinella_radians,(Globigerinella_siphonifera,(Globigerinella_calida,Globigerinella_adamsi))))),(((Trilobatus_sacculifer,Orbulina_universa),Sphaeroidinella_dehiscens),((Globoturborotalita_rubescens,Globoturborotalita_tenellus),((Globigerinoides_elongatus,Globigerinoides_conglobatus),(Globigerinoides_white,Globigerinoides_ruber)))))),(Hastigerina_pelagica,Hastigerinella_digitata)),(((((((Globorotalia_tumida,Globorotalia_ungulata),Globorotalia_menardii),(((Globorotalia_cavernula,Globorotalia_truncatulinoides),Globorotalia_crassaformis),(Globorotalia_scitula,(Globorotalia_theyeri,Globorotalia_hirsuta)))),Globoconella_inflata),((Pulleniatina_obliquiloculata,Neogloboquadrina_dutertrei),(Neogloboquadrina_incompta,Neogloboquadrina_pachyderma))),Globoquadrina_conglomerata),Globorotaloides_hexagonus,Berggrenia_pumilio,Dentigloborotalia_anfracta)),((((Candeina_nitida,Globigerinita_glutinata),Globigerinita_uvula),Globigerinita_minuta),(Tenuitella_fleisheri,Tenuitella_iota,Tenuitella_parkerae)));")
plankforams <- read.tree(text=";")

ggtree(plankforams, layout ='rectangular', ladderize=TRUE) + geom_tiplab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, color = "blue") +
  xlim(0,15)

write.nexus(plankforams, file = "data/phylogeny/kucera_tree.nex")
write.tree(plankforams, file = "data/phylogeny/kucera_tree.tre")

tree <- read.tree(file = "data/phylogeny/kucera_tree.tre")

tree <- compute.brlen(tree, 1) # sets all branch lengths equal to one
cophenetic(tree)
dist.nodes(tree)
branching.times(tree)
