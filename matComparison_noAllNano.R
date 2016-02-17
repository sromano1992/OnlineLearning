setwd('unisa/II_anno/Bioinformatica/progetto/myWork/')
library(grid)
library(gridExtra)
library(igraph)
attach(mtcars)
require(ggplot2)
library(fpc)
library(caret)
library(cluster)
con = gzcon(url('https://github.com/systematicinvestor/SIT/raw/master/sit.gz', 'rb'))
source(con)
close(con)
source("KTI.r")
source("plotTable.R")
source("MSC.R")
source("kendalltau.R")
source("../kendall_tao.R")
source("plotCommunities.R")
load("~/unisa/II_anno/Bioinformatica/progetto/myWork/workspace/.RData")

#ADJ = distance matrix without CLR 
#SubMatrix without some nano
nanoToRemoveIndex = sample(1:29, 10);
write("REMOVED", file = "removedNano.txt")
write(nanoToRemoveIndex, file = "removedNano.txt",ncolumns = length(nanoToRemoveIndex),append=TRUE)
write(nano[nanoToRemoveIndex], file = "removedNano.txt",ncolumns = length(nanoToRemoveIndex),append=TRUE)
write("USED", file = "removedNano.txt",append=TRUE)
nanoNewIndex = matrix(1:29);
nanoNewIndex = nanoNewIndex[!nanoNewIndex %in% nanoToRemoveIndex]
write(nanoNewIndex, file = "removedNano.txt",ncolumns = length(nanoNewIndex),append=TRUE)
write(nano[nanoNewIndex], file = "removedNano.txt",ncolumns = length(nanoNewIndex),append=TRUE)


newNanoNano = ADJ[nano,nano][nanoNewIndex,nanoNewIndex];
newNanoDrugs = ADJ[nano,drugs][nanoNewIndex,];
newNanoDisease = ADJ[nano,disease][nanoNewIndex,];
newNanoChemical = ADJ[nano,chemical][nanoNewIndex,];
ADJ_no_all_nano = rbind(cbind(newNanoNano,
                              newNanoDrugs,
                              newNanoDisease,
                              newNanoChemical),
                        cbind(t(newNanoDrugs),
                              ADJ[drugs,drugs],
                              ADJ[drugs,disease],
                              ADJ[drugs,chemical]),
                        cbind(t(newNanoDisease),
                              ADJ[disease,drugs],
                              ADJ[disease,disease],
                              ADJ[disease,chemical]),
                        cbind(t(newNanoChemical),
                              ADJ[chemical,drugs],
                              ADJ[chemical,disease],
                              ADJ[chemical,chemical]));
#dump("KTI", file='KTI.r');
#Computing CLR
#W_ADJ = KTI(KTDM_mat = abs(ADJ),SIGN = sign(ADJ))
W_ADJ_no_all_nano = KTI(KTDM_mat = abs(ADJ_no_all_nano),SIGN = sign(ADJ_no_all_nano))

#Getting submatrices to compare
W_ADJ_nano_nano = W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]];
W_ADJ_nano_drugs = W_ADJ[nano[nanoNewIndex],drugs];
W_ADJ_nano_disease = W_ADJ[nano[nanoNewIndex],disease];
W_ADJ_nano_chemical = W_ADJ[nano[nanoNewIndex],chemical];
W_ADJ_drugs_nano = W_ADJ[drugs,nano[nanoNewIndex]];
W_ADJ_drugs_drugs = W_ADJ[drugs,drugs];
W_ADJ_drugs_disease = W_ADJ[drugs,disease];
W_ADJ_drugs_chemical = W_ADJ[drugs,chemical];
W_ADJ_disease_nano = W_ADJ[disease, nano[nanoNewIndex]];
W_ADJ_disease_drugs = W_ADJ[disease, drugs];
W_ADJ_disease_disease = W_ADJ[disease, disease];
W_ADJ_disease_chemical = W_ADJ[disease, chemical];
W_ADJ_chemical_nano = W_ADJ[chemical, nano[nanoNewIndex]];
W_ADJ_chemical_drugs = W_ADJ[chemical, drugs];
W_ADJ_chemical_disease = W_ADJ[chemical, disease];
W_ADJ_chemical_chemical = W_ADJ[chemical, chemical];
W_ADJ_nanoDrugs_nanoDrugs = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],drugs]),
  cbind(W_ADJ[drugs,nano[nanoNewIndex]],W_ADJ[drugs,drugs])
);
W_ADJ_nanoDisease_nanoDisease = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],disease]),
  cbind(W_ADJ[disease,nano[nanoNewIndex]],W_ADJ[disease,disease])
);
W_ADJ_nanoChemical_nanoChemical = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ[chemical,nano[nanoNewIndex]],W_ADJ[chemical,chemical])
);
W_ADJ_drugsDisease_drugsDisease = rbind(
  cbind(W_ADJ[drugs,drugs],W_ADJ[drugs,disease]),
  cbind(W_ADJ[disease,drugs],W_ADJ[disease,disease])
);
W_ADJ_drugsChemical_drugsChemical = rbind(
  cbind(W_ADJ[drugs,drugs],W_ADJ[drugs,chemical]),
  cbind(W_ADJ[chemical,drugs],W_ADJ[chemical,chemical])
);
W_ADJ_diseaseChemical_diseaseChemical = rbind(
  cbind(W_ADJ[disease,disease],W_ADJ[disease,chemical]),
  cbind(W_ADJ[chemical,disease],W_ADJ[chemical,chemical])
);
W_ADJ_nanoDrugsDisease = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],drugs],W_ADJ[nano[nanoNewIndex],disease]),
  cbind(W_ADJ[drugs,nano[nanoNewIndex]],W_ADJ[drugs,drugs],W_ADJ[drugs,disease]),
  cbind(W_ADJ[disease,nano[nanoNewIndex]],W_ADJ[disease,drugs],W_ADJ[disease,disease])
);
W_ADJ_nanoDiseaseChemical = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],disease],W_ADJ[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ[disease,nano[nanoNewIndex]],W_ADJ[disease,disease],W_ADJ[disease,chemical]),
  cbind(W_ADJ[chemical,nano[nanoNewIndex]],W_ADJ[chemical,disease],W_ADJ[chemical,chemical])
);
W_ADJ_drugsDiseaseChemical = rbind(
  cbind(W_ADJ[drugs,drugs],W_ADJ[drugs,disease],W_ADJ[drugs,chemical]),
  cbind(W_ADJ[disease,drugs],W_ADJ[disease,disease],W_ADJ[disease,chemical]),
  cbind(W_ADJ[chemical,drugs],W_ADJ[chemical,disease],W_ADJ[chemical,chemical])
);
W_ADJ_nanoDrugsChemical = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],drugs],W_ADJ[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ[drugs,nano[nanoNewIndex]],W_ADJ[drugs,drugs],W_ADJ[drugs,chemical]),
  cbind(W_ADJ[chemical,nano[nanoNewIndex]],W_ADJ[chemical,drugs],W_ADJ[chemical,chemical])
);
W_ADJ_sub = rbind(
  cbind(W_ADJ[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ[nano[nanoNewIndex],drugs],W_ADJ[nano[nanoNewIndex],disease],W_ADJ[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ[drugs,nano[nanoNewIndex]],W_ADJ[drugs,drugs],W_ADJ[drugs,disease],W_ADJ[drugs,chemical]),
  cbind(W_ADJ[disease,nano[nanoNewIndex]],W_ADJ[disease,drugs],W_ADJ[disease,disease],W_ADJ[disease,chemical]),
  cbind(W_ADJ[chemical,nano[nanoNewIndex]],W_ADJ[chemical,drugs],W_ADJ[chemical,disease],W_ADJ[chemical,chemical])
);

W_ADJ_no_all_nano_nano_nano = W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]];
W_ADJ_no_all_nano_nano_drugs = W_ADJ_no_all_nano[nano[nanoNewIndex],drugs];
W_ADJ_no_all_nano_nano_disease = W_ADJ_no_all_nano[nano[nanoNewIndex],disease];
W_ADJ_no_all_nano_nano_chemical = W_ADJ_no_all_nano[nano[nanoNewIndex],chemical];
W_ADJ_no_all_nano_drugs_nano = W_ADJ_no_all_nano[drugs,nano[nanoNewIndex]];
W_ADJ_no_all_nano_drugs_drugs = W_ADJ_no_all_nano[drugs,drugs];
W_ADJ_no_all_nano_drugs_disease = W_ADJ_no_all_nano[drugs,disease];
W_ADJ_no_all_nano_drugs_chemical = W_ADJ_no_all_nano[drugs,chemical];
W_ADJ_no_all_nano_disease_nano = W_ADJ_no_all_nano[disease, nano[nanoNewIndex]];
W_ADJ_no_all_nano_disease_drugs = W_ADJ_no_all_nano[disease, drugs];
W_ADJ_no_all_nano_disease_disease = W_ADJ_no_all_nano[disease, disease];
W_ADJ_no_all_nano_disease_chemical = W_ADJ_no_all_nano[disease, chemical];
W_ADJ_no_all_nano_chemical_chemical = W_ADJ_no_all_nano[chemical, chemical];
W_ADJ_no_all_nano_nanoDrugs_nanoDrugs = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],drugs]),
  cbind(W_ADJ_no_all_nano[drugs,nano[nanoNewIndex]],W_ADJ_no_all_nano[drugs,drugs])
);
W_ADJ_no_all_nano_nanoDisease_nanoDisease = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],disease]),
  cbind(W_ADJ_no_all_nano[disease,nano[nanoNewIndex]],W_ADJ_no_all_nano[disease,disease])
);
W_ADJ_no_all_nano_nanoChemical_nanoChemical = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ_no_all_nano[chemical,nano[nanoNewIndex]],W_ADJ_no_all_nano[chemical,chemical])
);
W_ADJ_no_all_nano_drugsDisease_drugsDisease = rbind(
  cbind(W_ADJ_no_all_nano[drugs,drugs],W_ADJ_no_all_nano[drugs,disease]),
  cbind(W_ADJ_no_all_nano[disease,drugs],W_ADJ_no_all_nano[disease,disease])
);
W_ADJ_no_all_nano_drugsChemical_drugsChemical = rbind(
  cbind(W_ADJ_no_all_nano[drugs,drugs],W_ADJ_no_all_nano[drugs,chemical]),
  cbind(W_ADJ_no_all_nano[chemical,drugs],W_ADJ_no_all_nano[chemical,chemical])
);
W_ADJ_no_all_nano_diseaseChemical_diseaseChemical = rbind(
  cbind(W_ADJ_no_all_nano[disease,disease],W_ADJ_no_all_nano[disease,chemical]),
  cbind(W_ADJ_no_all_nano[chemical,disease],W_ADJ_no_all_nano[chemical,chemical])
);
W_ADJ_no_all_nano_nanoDrugsDisease = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],drugs],W_ADJ_no_all_nano[nano[nanoNewIndex],disease]),
  cbind(W_ADJ_no_all_nano[drugs,nano[nanoNewIndex]],W_ADJ_no_all_nano[drugs,drugs],W_ADJ_no_all_nano[drugs,disease]),
  cbind(W_ADJ_no_all_nano[disease,nano[nanoNewIndex]],W_ADJ_no_all_nano[disease,drugs],W_ADJ_no_all_nano[disease,disease])
);
W_ADJ_no_all_nano_nanoDiseaseChemical = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],disease],W_ADJ_no_all_nano[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ_no_all_nano[disease,nano[nanoNewIndex]],W_ADJ_no_all_nano[disease,disease],W_ADJ_no_all_nano[disease,chemical]),
  cbind(W_ADJ_no_all_nano[chemical,nano[nanoNewIndex]],W_ADJ_no_all_nano[chemical,disease],W_ADJ_no_all_nano[chemical,chemical])
);
W_ADJ_no_all_nano_drugsDiseaseChemical = rbind(
  cbind(W_ADJ_no_all_nano[drugs,drugs],W_ADJ_no_all_nano[drugs,disease],W_ADJ_no_all_nano[drugs,chemical]),
  cbind(W_ADJ_no_all_nano[disease,drugs],W_ADJ_no_all_nano[disease,disease],W_ADJ_no_all_nano[disease,chemical]),
  cbind(W_ADJ_no_all_nano[chemical,drugs],W_ADJ_no_all_nano[chemical,disease],W_ADJ_no_all_nano[chemical,chemical])
);
W_ADJ_no_all_nano_nanoDrugsChemical = rbind(
  cbind(W_ADJ_no_all_nano[nano[nanoNewIndex],nano[nanoNewIndex]],W_ADJ_no_all_nano[nano[nanoNewIndex],drugs],W_ADJ_no_all_nano[nano[nanoNewIndex],chemical]),
  cbind(W_ADJ_no_all_nano[drugs,nano[nanoNewIndex]],W_ADJ_no_all_nano[drugs,drugs],W_ADJ_no_all_nano[drugs,chemical]),
  cbind(W_ADJ_no_all_nano[chemical,nano[nanoNewIndex]],W_ADJ_no_all_nano[chemical,drugs],W_ADJ_no_all_nano[chemical,chemical])
);

# #Nano-Nano analysis plot
# tmp = W_ADJ_nano_nano;
# diag(tmp) = 0;
# g = graph.adjacency(tmp,mode="undirected",weighted = T)
# plot(g,edge.width = E(g)$weight * 1.2)
# 
# dev.new()
# tmp = W_ADJ_no_all_nano_nano_nano;
# diag(tmp) = 0;
# g1 = graph.adjacency(tmp,mode="undirected",weighted = T)
# plot(g1,edge.width = E(g1)$weight * 1.2)
# 
# par(mfrow=c(1,2))
# plot(g,edge.width = E(g)$weight * 1.2, main="Nano_nano complete CLR")
# plot(g1,edge.width = E(g1)$weight * 1.2, main="Nano_nano CLR without nano")

#Cluster analysis
k = 4;
iterations = 1000000;
#find k
#diss <- as.matrix(dist(W_ADJ_nano_nano));
#for (i in 4:28){
#  cluster_W_ADJ_nano_nano <- kmeans(W_ADJ_nano_nano,i,iterations)
#  sil <-mean(silhouette(as.integer(cluster_W_ADJ_nano_nano$cluster),dmatrix=diss)[,3])
#  print(sil)
#}
#nano-nano
cluster_W_ADJ_nano_nano <- kmeans(W_ADJ_nano_nano,k,iterations)
cluster_W_ADJ_no_all_nano_nano_nano <- kmeans(W_ADJ_no_all_nano_nano_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_nano$cluster,cluster_W_ADJ_no_all_nano_nano_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot <- ggplot(confusion)
plot + ggtitle("Nano-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-drugs
cluster_W_ADJ_nano_drugs <- kmeans(W_ADJ_nano_drugs,k,iterations)
cluster_W_ADJ_no_all_nano_nano_drugs <- kmeans(W_ADJ_no_all_nano_nano_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_drugs$cluster,cluster_W_ADJ_no_all_nano_nano_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Nano-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-disease
cluster_W_ADJ_nano_disease <- kmeans(W_ADJ_nano_disease,k,iterations)
cluster_W_ADJ_no_all_nano_nano_disease <- kmeans(W_ADJ_no_all_nano_nano_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_disease$cluster,cluster_W_ADJ_no_all_nano_nano_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Nano-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-chemical
cluster_W_ADJ_nano_chemical <- kmeans(W_ADJ_nano_chemical,k,iterations)
cluster_W_ADJ_no_all_nano_nano_chemical <- kmeans(W_ADJ_no_all_nano_nano_chemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_chemical$cluster,cluster_W_ADJ_no_all_nano_nano_chemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Nano-Chemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-nano
cluster_W_ADJ_drugs_nano <- kmeans(W_ADJ_drugs_nano,k,iterations)
cluster_W_ADJ_no_all_nano_drugs_nano <- kmeans(W_ADJ_no_all_nano_drugs_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_nano$cluster,cluster_W_ADJ_no_all_nano_drugs_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-drugs
cluster_W_ADJ_drugs_drugs <- kmeans(W_ADJ_drugs_drugs,k,iterations)
cluster_W_ADJ_no_all_nano_drugs_drugs <- kmeans(W_ADJ_no_all_nano_drugs_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_drugs$cluster,cluster_W_ADJ_no_all_nano_drugs_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-disease
cluster_W_ADJ_drugs_disease <- kmeans(W_ADJ_drugs_disease,k,iterations)
cluster_W_ADJ_no_all_nano_drugs_disease <- kmeans(W_ADJ_no_all_nano_drugs_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_disease$cluster,cluster_W_ADJ_no_all_nano_drugs_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-chemical
cluster_W_ADJ_drugs_chemical <- kmeans(W_ADJ_drugs_chemical,k,iterations)
cluster_W_ADJ_no_all_nano_drugs_chemical <- kmeans(W_ADJ_no_all_nano_drugs_chemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_chemical$cluster,cluster_W_ADJ_no_all_nano_drugs_chemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Chemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-nano
cluster_W_ADJ_disease_nano <- kmeans(W_ADJ_disease_nano,k,iterations)
cluster_W_ADJ_no_all_nano_disease_nano <- kmeans(W_ADJ_no_all_nano_disease_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_nano$cluster,cluster_W_ADJ_no_all_nano_disease_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-drugs
cluster_W_ADJ_disease_drugs <- kmeans(W_ADJ_disease_drugs,k,iterations)
cluster_W_ADJ_no_all_nano_disease_drugs <- kmeans(W_ADJ_no_all_nano_disease_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_drugs$cluster,cluster_W_ADJ_no_all_nano_disease_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-disease
cluster_W_ADJ_disease_disease <- kmeans(W_ADJ_disease_disease,k,iterations)
cluster_W_ADJ_no_all_nano_disease_disease <- kmeans(W_ADJ_no_all_nano_disease_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_disease$cluster,cluster_W_ADJ_no_all_nano_disease_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-chemical
cluster_W_ADJ_disease_chemical <- kmeans(W_ADJ_disease_chemical,k,iterations)
cluster_W_ADJ_no_all_nano_disease_chemical <- kmeans(W_ADJ_no_all_nano_disease_chemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_chemical$cluster,cluster_W_ADJ_no_all_nano_disease_chemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Chemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#chemical-chemical
cluster_W_ADJ_chemical_chemical <- kmeans(W_ADJ_chemical_chemical,k,iterations)
cluster_W_ADJ_no_all_nano_chemical_chemical <- kmeans(W_ADJ_no_all_nano_chemical_chemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_chemical_chemical$cluster,cluster_W_ADJ_no_all_nano_chemical_chemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Chemical-Chemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nanoDrugs-nanoDrugs
cluster_W_ADJ_nanoDrugs_nanoDrugs <- kmeans(W_ADJ_nanoDrugs_nanoDrugs,k,iterations)
cluster_W_ADJ_no_all_nano_nanoDrugs_nanoDrugs <- kmeans(W_ADJ_no_all_nano_nanoDrugs_nanoDrugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDrugs_nanoDrugs$cluster,cluster_W_ADJ_no_all_nano_nanoDrugs_nanoDrugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDrugs-NanoDrugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nanoDisease-nanoDisease
cluster_W_ADJ_nanoDisease_nanoDisease <- kmeans(W_ADJ_nanoDisease_nanoDisease,k,iterations)
cluster_W_ADJ_no_all_nano_nanoDisease_nanoDisease <- kmeans(W_ADJ_no_all_nano_nanoDisease_nanoDisease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDisease_nanoDisease$cluster,cluster_W_ADJ_no_all_nano_nanoDisease_nanoDisease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDisease-NanoDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nanoChemical-nanoChemical
cluster_W_ADJ_nanoChemical_nanoChemical <- kmeans(W_ADJ_nanoChemical_nanoChemical,k,iterations)
cluster_W_ADJ_no_all_nano_nanoChemical_nanoChemical <- kmeans(W_ADJ_no_all_nano_nanoChemical_nanoChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoChemical_nanoChemical$cluster,cluster_W_ADJ_no_all_nano_nanoChemical_nanoChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoChemical-NanoChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugsDisease-drugsDisease
cluster_W_ADJ_drugsDisease_drugsDisease <- kmeans(W_ADJ_drugsDisease_drugsDisease,k,iterations)
cluster_W_ADJ_no_all_nano_drugsDisease_drugsDisease <- kmeans(W_ADJ_no_all_nano_drugsDisease_drugsDisease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugsDisease_drugsDisease$cluster,cluster_W_ADJ_no_all_nano_drugsDisease_drugsDisease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("DrugsDisease-DrugsDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugsChemical-drugsChemical
cluster_W_ADJ_drugsChemical_drugsChemical <- kmeans(W_ADJ_drugsChemical_drugsChemical,k,iterations)
cluster_W_ADJ_no_all_nano_drugsChemical_drugsChemical <- kmeans(W_ADJ_no_all_nano_drugsChemical_drugsChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugsChemical_drugsChemical$cluster,cluster_W_ADJ_no_all_nano_drugsChemical_drugsChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("DrugsChemical-DrugsChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#diseaseChemical-diseaseChemical
cluster_W_ADJ_diseaseChemical_diseaseChemical <- kmeans(W_ADJ_diseaseChemical_diseaseChemical,k,iterations)
cluster_W_ADJ_no_all_nano_diseaseChemical_diseaseChemical <- kmeans(W_ADJ_no_all_nano_diseaseChemical_diseaseChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_diseaseChemical_diseaseChemical$cluster,cluster_W_ADJ_no_all_nano_diseaseChemical_diseaseChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("DiseaseChemical-DiseaseChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-drugs-disease
cluster_W_ADJ_nanoDrugsDisease <- kmeans(W_ADJ_nanoDrugsDisease,k,iterations)
cluster_W_ADJ_no_all_nano <- kmeans(W_ADJ_no_all_nano_nanoDrugsDisease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDrugsDisease$cluster,cluster_W_ADJ_no_all_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDrugsDisease-NanoDrugsDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-disease-chemical
cluster_W_ADJ_nanoDiseaseChemical <- kmeans(W_ADJ_nanoDiseaseChemical,k,iterations)
cluster_W_ADJ_no_all_nano_nanoDiseaseChemical <- kmeans(W_ADJ_no_all_nano_nanoDiseaseChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDiseaseChemical$cluster,cluster_W_ADJ_no_all_nano_nanoDiseaseChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDiseaseChemical-NanoDiseaseChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-disease-chemical
cluster_W_ADJ_drugsDiseaseChemical <- kmeans(W_ADJ_drugsDiseaseChemical,k,iterations)
cluster_W_ADJ_no_all_nano_drugsDiseaseChemical <- kmeans(W_ADJ_no_all_nano_drugsDiseaseChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugsDiseaseChemical$cluster,cluster_W_ADJ_no_all_nano_drugsDiseaseChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("DrugsDiseaseChemical-DrugsDiseaseChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-drugs-chemical
cluster_W_ADJ_nanoDugsChemical <- kmeans(W_ADJ_nanoDrugsChemical,k,iterations)
cluster_W_ADJ_no_all_nano_nanoDrugsChemical <- kmeans(W_ADJ_no_all_nano_nanoDrugsChemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDugsChemical$cluster,cluster_W_ADJ_no_all_nano_nanoDrugsChemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDrugsChemical-NanoDrugsChemical clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#All
cluster_W_ADJ <- kmeans(W_ADJ_sub,k,iterations)
cluster_W_ADJ_no_all_nano <- kmeans(W_ADJ_no_all_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ$cluster,cluster_W_ADJ_no_all_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("All data clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 

#Social network analysis test
#Nano-Nano analysis
g_nano_nano = graph.adjacency(W_ADJ_nano_nano,mode="undirected",weighted = T)
g_no_all_nano_nano_nano = graph.adjacency(W_ADJ_no_all_nano_nano_nano,mode="undirected",weighted = T)
degree_nano_nano = degree(g_nano_nano);
degree_no_all_nano_nano_nano = degree(g_no_all_nano_nano_nano);
tmp=sort(degree_nano_nano,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nano_nano,decreasing = TRUE);
degreeNanoNano_distance = kendallTauDistance(names(tmp), names(tmp1));
table = rbind(cbind('KendallTau Distance',degreeNanoNano_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,'degree_nanoNano.pdf');
# plot(degree_nano_nano, degree_no_all_nano_nano_nano);
# plot(g_nano_nano,
#      vertex.size = round(degree_nano_nano, 2),
#      vertex.color = rainbow(vcount(g_nano_nano)))
# dev.new()
# plot(g_no_all_nano_nano_nano,
#      vertex.size = round(degree_no_all_nano_nano_nano, 2),
#      vertex.color = rainbow(vcount(g_no_all_nano_nano_nano)))
#Drugs-Drugs analysis
g_drugs_drugs = graph.adjacency(W_ADJ_drugs_drugs,mode="undirected",weighted = T)
g_no_all_nano_drugs_drugs = graph.adjacency(W_ADJ_no_all_nano_drugs_drugs,mode="undirected",weighted = T)
degree_drugs_drugs = degree(g_drugs_drugs);
degree_no_all_nano_drugs_drugs = degree(g_no_all_nano_drugs_drugs);
tmp=sort(degree_drugs_drugs,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_drugs_drugs,decreasing = TRUE);
degreeDrugsDrugs_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsDrugs_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_drugsDrugs.pdf");
#Disease-Disease analysis
g_disease_disease = graph.adjacency(W_ADJ_disease_disease,mode="undirected",weighted = T)
g_no_all_nano_disease_disease = graph.adjacency(W_ADJ_no_all_nano_disease_disease,mode="undirected",weighted = T)
degree_disease_disease = degree(g_disease_disease);
degree_no_all_nano_disease_disease = degree(g_no_all_nano_disease_disease);
tmp=sort(degree_disease_disease,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_disease_disease,decreasing = TRUE);
degreeDiseaseDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDiseaseDisease_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_diseaseDisease.pdf");
#Chemical-Chemical analysis
g_chemical_chemical = graph.adjacency(W_ADJ_chemical_chemical,mode="undirected",weighted = T)
g_no_all_nano_chemical_chemical = graph.adjacency(W_ADJ_no_all_nano_chemical_chemical,mode="undirected",weighted = T)
degree_chemical_chemical = degree(g_chemical_chemical);
degree_no_all_nano_chemical_chemical = degree(g_no_all_nano_chemical_chemical);
tmp=sort(degree_chemical_chemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_chemical_chemical,decreasing = TRUE);
degreeChemicalChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeChemicalChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_chemicalChemical.pdf");
#NanoDrugs-NanoDrugs analysis
g_nanoDrugs_nanoDrugs = graph.adjacency(W_ADJ_nanoDrugs_nanoDrugs,mode="undirected",weighted = T)
g_no_all_nano_nanoDrugs_nanoDrugs = graph.adjacency(W_ADJ_no_all_nano_nanoDrugs_nanoDrugs,mode="undirected",weighted = T)
degree_nanoDrugs_nanoDrugs = degree(g_nanoDrugs_nanoDrugs);
degree_no_all_nano_nanoDrugs_nanoDrugs = degree(g_no_all_nano_nanoDrugs_nanoDrugs);
tmp=sort(degree_nanoDrugs_nanoDrugs,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nanoDrugs_nanoDrugs,decreasing = TRUE);
degreeNanoDrugs_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDrugs_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_nanoDrugs_nanoDrugs.pdf");
#NanoDisease-NanoDisease analysis
g_nanoDisease_nanoDisease = graph.adjacency(W_ADJ_nanoDisease_nanoDisease,mode="undirected",weighted = T)
g_no_all_nano_nanoDisease_nanoDisease = graph.adjacency(W_ADJ_no_all_nano_nanoDisease_nanoDisease,mode="undirected",weighted = T)
degree_nanoDisease_nanoDisease = degree(g_nanoDisease_nanoDisease);
degree_no_all_nano_nanoDisease_nanoDisease = degree(g_no_all_nano_nanoDisease_nanoDisease);
tmp=sort(degree_nanoDisease_nanoDisease,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nanoDisease_nanoDisease,decreasing = TRUE);
degreeNanoDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDisease_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_nanoDisease_nanoDisease.pdf");
#NanoChemical-NanoChemical analysis
g_nanoChemical_nanoChemical = graph.adjacency(W_ADJ_nanoChemical_nanoChemical,mode="undirected",weighted = T)
g_no_all_nano_nanoChemical_nanoChemical = graph.adjacency(W_ADJ_no_all_nano_nanoChemical_nanoChemical,mode="undirected",weighted = T)
degree_nanoChemical_nanoChemical = degree(g_nanoChemical_nanoChemical);
degree_no_all_nano_nanoChemical_nanoChemical = degree(g_no_all_nano_nanoChemical_nanoChemical);
tmp=sort(degree_nanoChemical_nanoChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nanoChemical_nanoChemical,decreasing = TRUE);
degreeNanoChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_nanoChemical_nanoChemical.pdf");
#DrugDisease-DrugDisease analysis
g_drugsDisease_drugsDisease = graph.adjacency(W_ADJ_drugsDisease_drugsDisease,mode="undirected",weighted = T)
g_no_all_nano_drugsDisease_drugsDisease = graph.adjacency(W_ADJ_no_all_nano_drugsDisease_drugsDisease,mode="undirected",weighted = T)
degree_drugsDisease_drugsDisease = degree(g_drugsDisease_drugsDisease);
degree_no_all_nano_drugsDisease_drugsDisease = degree(g_no_all_nano_drugsDisease_drugsDisease);
tmp=sort(degree_drugsDisease_drugsDisease,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_drugsDisease_drugsDisease,decreasing = TRUE);
degreeDrugsDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsDisease_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_drugsDisease_drugsDisease.pdf");
#DrugChemical-DrugChemical analysis
g_drugsChemical_drugsChemical = graph.adjacency(W_ADJ_drugsChemical_drugsChemical,mode="undirected",weighted = T)
g_no_all_nano_drugsChemical_drugsChemical = graph.adjacency(W_ADJ_no_all_nano_drugsChemical_drugsChemical,mode="undirected",weighted = T)
degree_drugsChemical_drugsChemical = degree(g_drugsChemical_drugsChemical);
degree_no_all_nano_drugsChemical_drugsChemical = degree(g_no_all_nano_drugsChemical_drugsChemical);
tmp=sort(degree_drugsChemical_drugsChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_drugsChemical_drugsChemical,decreasing = TRUE);
degreeDrugsChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_drugsChemical_drugsChemical.pdf");
#DiseaseChemical-DiseaseChemical analysis
g_diseaseChemical_diseaseChemical = graph.adjacency(W_ADJ_diseaseChemical_diseaseChemical,mode="undirected",weighted = T)
g_no_all_nano_diseaseChemical_diseaseChemical = graph.adjacency(W_ADJ_no_all_nano_diseaseChemical_diseaseChemical,mode="undirected",weighted = T)
degree_diseaseChemical_diseaseChemical = degree(g_diseaseChemical_diseaseChemical);
degree_no_all_nano_diseaseChemical_diseaseChemical = degree(g_no_all_nano_diseaseChemical_diseaseChemical);
tmp=sort(degree_diseaseChemical_diseaseChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_diseaseChemical_diseaseChemical,decreasing = TRUE);
degreeDiseaseChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDiseaseChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_diseaseChemical_diseaseChemical.pdf");
#nanoDrugsDisease
g_nanoDrugsDisease_nanoDrugsDisease = graph.adjacency(W_ADJ_nanoDrugsDisease,mode="undirected",weighted = T)
g_no_all_nano_nanoDrugsDisease_nanoDrugsDisease = graph.adjacency(W_ADJ_no_all_nano_nanoDrugsDisease,mode="undirected",weighted = T)
degree_nanoDrugsDisease_nanoDrugsDisease = degree(g_nanoDrugsDisease_nanoDrugsDisease);
degree_no_all_nano_diseaseChemical_diseaseChemical = degree(g_no_all_nano_nanoDrugsDisease_nanoDrugsDisease);
tmp=sort(degree_nanoDrugsDisease_nanoDrugsDisease,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_diseaseChemical_diseaseChemical,decreasing = TRUE);
degreeNanoDrugsDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDrugsDisease_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_nanoDrugsDisease_nanoDrugsDisease.pdf");
#nanoDiseaseChamical
g_nanoDiseaseChemical_nanoDiseaseChemical = graph.adjacency(W_ADJ_nanoDiseaseChemical,mode="undirected",weighted = T)
g_no_all_nano_nanoDiseaseChemical_nanoDiseaseChemical = graph.adjacency(W_ADJ_no_all_nano_nanoDiseaseChemical,mode="undirected",weighted = T)
degree_nanoDiseaseChemical_nanoDiseaseChemical = degree(g_nanoDiseaseChemical_nanoDiseaseChemical);
degree_no_all_nano_nanoDiseaseChemical_nanoDiseaseChemical = degree(g_no_all_nano_nanoDiseaseChemical_nanoDiseaseChemical);
tmp=sort(degree_nanoDiseaseChemical_nanoDiseaseChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nanoDiseaseChemical_nanoDiseaseChemical,decreasing = TRUE);
degreeNanoDiseaseChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDiseaseChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_nanoDiseaseChemical_nanoDiseaseChemical.pdf");
#drugsDiseaseChemical
g_drugsDiseaseChemical_drugsDiseaseChemical = graph.adjacency(W_ADJ_drugsDiseaseChemical,mode="undirected",weighted = T)
g_no_all_nano_drugsDiseaseChemical_drugsDiseaseChemical = graph.adjacency(W_ADJ_no_all_nano_drugsDiseaseChemical,mode="undirected",weighted = T)
degree_drugsDiseaseChemical_drugsDiseaseChemical = degree(g_drugsDiseaseChemical_drugsDiseaseChemical);
degree_no_all_nano_drugsDiseaseChemical_drugsDiseaseChemical = degree(g_no_all_nano_drugsDiseaseChemical_drugsDiseaseChemical);
tmp=sort(degree_drugsDiseaseChemical_drugsDiseaseChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_drugsDiseaseChemical_drugsDiseaseChemical,decreasing = TRUE);
degreeDrugsDiseaseChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsDiseaseChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_drugsDiseaseChemical_drugsDiseaseChemical.pdf");
#nanoDrugsChemical
g_nanoDrugsChemical_nanoDrugsChemical = graph.adjacency(W_ADJ_nanoDrugsChemical,mode="undirected",weighted = T)
g_no_all_nano_nanoDrugsChemical_nanoDrugsChemical = graph.adjacency(W_ADJ_no_all_nano_nanoDrugsChemical,mode="undirected",weighted = T)
degree_nanoDrugsChemical_nanoDrugsChemical = degree(g_nanoDrugsChemical_nanoDrugsChemical);
degree_no_all_nano_nanoDrugsChemical_nanoDrugsChemical = degree(g_no_all_nano_nanoDrugsChemical_nanoDrugsChemical);
tmp=sort(degree_nanoDrugsChemical_nanoDrugsChemical,decreasing = TRUE);
tmp1=sort(degree_no_all_nano_nanoDrugsChemical_nanoDrugsChemical,decreasing = TRUE);
degreeNanoDrugsChemical_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDrugsChemical_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_nanoDrugsChemical_nanoDrugsChemical.pdf");
#all matrix
g_W_ADJ_no_all_nano = graph.adjacency(W_ADJ_no_all_nano,mode="undirected",weighted = T)
g_W_ADJ_sub = graph.adjacency(W_ADJ_sub,mode="undirected",weighted = T)
degree_W_ADJ_sub = degree(g_W_ADJ_sub);
degree_W_ADJ_no_all_nano = degree(g_W_ADJ_no_all_nano);
tmp=sort(degree_W_ADJ_sub,decreasing = TRUE);
tmp1=sort(degree_W_ADJ_no_all_nano,decreasing = TRUE);
degreeAll_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeAll_distance,'',''),cbind('DEGREE (no nano)','NAME','DEGREE (with nano)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
table[,] = substring(table[,],1,20)
plotTable(table,"degree_All.pdf");

#Communities check
#nano-nano
# par(mfrow=c(1,2))
# plot(walktrap.community(g_nano_nano),simplify(g_nano_nano,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_nano_nano),simplify(g_no_all_nano_nano_nano,remove.loops = TRUE))
# #drugs-drugs
# par(mfrow=c(1,2))
# plot(walktrap.community(g_drugs_drugs),simplify(g_drugs_drugs,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_drugs_drugs),simplify(g_no_all_nano_drugs_drugs,remove.loops = TRUE))
# #disease-disease
# par(mfrow=c(1,2))
# plot(walktrap.community(g_disease_disease),simplify(g_disease_disease,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_disease_disease),simplify(g_no_all_nano_disease_disease,remove.loops = TRUE))
# #nanoDrugs-nanoDrugs
# par(mfrow=c(1,2))
# plot(walktrap.community(g_nanoDrugs_nanoDrugs),simplify(g_nanoDrugs_nanoDrugs,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_nanoDrugs_nanoDrugs),simplify(g_no_all_nano_nanoDrugs_nanoDrugs,remove.loops = TRUE))
# #nanoDisease-nanoDisease
# par(mfrow=c(1,2))
# plot(walktrap.community(g_nanoDisease_nanoDisease),simplify(g_nanoDisease_nanoDisease,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_nanoDisease_nanoDisease),simplify(g_no_all_nano_nanoDisease_nanoDisease,remove.loops = TRUE))
# #drugsDisease-drugsDisease
# par(mfrow=c(1,2))
# plot(walktrap.community(g_drugsDisease_drugsDisease),simplify(g_drugsDisease_drugsDisease,remove.loops = TRUE))
# plot(walktrap.community(g_no_all_nano_drugsDisease_drugsDisease),simplify(g_no_all_nano_drugsDisease_drugsDisease,remove.loops = TRUE))
# #nano-drugs-disease
# par(mfrow=c(1,2))
# plot(walktrap.community(g_W_ADJ_sub),simplify(g_W_ADJ_sub,remove.loops = TRUE))
# plot(walktrap.community(g_W_ADJ_no_all_nano),simplify(g_W_ADJ_no_all_nano,remove.loops = TRUE))

#Interactive plot
#nano-nano
tmp = g_nano_nano;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
edgeColor = matrix(0,dim(W_ADJ_nano_nano)[1],dim(W_ADJ_nano_nano)[2]);
edgeColor[W_ADJ_nano_nano<0]='red';
edgeColor[W_ADJ_nano_nano>=0]='green';
#tkplot(simplify(g_nano_nano,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)], edge.color=edgeColor,edge.curved =TRUE)
#gd <- get.data.frame(g_nano_nano, what = "edges")
#simpleNetwork(gd, fontSize = 12)
tmp = g_no_all_nano_nano_nano;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
edgeColor = matrix(0,dim(W_ADJ_no_all_nano_nano_nano)[1],dim(W_ADJ_no_all_nano_nano_nano)[2]);
edgeColor[W_ADJ_no_all_nano_nano_nano<0]='red';
edgeColor[W_ADJ_no_all_nano_nano_nano>=0]='green';
#tkplot(simplify(g_no_all_nano_nano_nano,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)], edge.color=edgeColor)
nano_nano_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoNano",nano_nano_communityDistance);
#drugs-drugs
tmp = g_drugs_drugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_drugs_drugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_all_nano_drugs_drugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_all_nano_drugs_drugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
drugs_drugs_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DrugsDrugs",drugs_drugs_communityDistance);
#disease-disease
tmp = g_disease_disease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_disease_disease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_all_nano_disease_disease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_all_nano_disease_disease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
disease_disease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DiseaseDisease",disease_disease_communityDistance);
#chemical-chemical
tmp = g_chemical_chemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
tmp = g_no_all_nano_chemical_chemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
chemical_chemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_ChemicalChemical",chemical_chemical_communityDistance);
#nanoDrugs-nanoDrugs
tmp = g_nanoDrugs_nanoDrugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_nanoDrugs_nanoDrugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_all_nano_nanoDrugs_nanoDrugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_all_nano_nanoDrugs_nanoDrugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
nanoDrugs_nanoDrugs_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDrugs",nanoDrugs_nanoDrugs_communityDistance);
#nanoDisease-nanoDisease
tmp = g_nanoDisease_nanoDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_nanoDisease_nanoDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_all_nano_nanoDisease_nanoDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_all_nano_nanoDisease_nanoDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
nanoDisease_nanoDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDisease",nanoDisease_nanoDisease_communityDistance);
#nanoChemical-nanoChemical
tmp = g_nanoChemical_nanoChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
tmp = g_no_all_nano_nanoChemical_nanoChemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
nanoChemical_nanoChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoChemical",nanoChemical_nanoChemical_communityDistance);
#drugsDisease-drugsDisease
tmp = g_drugsDisease_drugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_drugsDisease_drugsDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_all_nano_drugsDisease_drugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_all_nano_drugsDisease_drugsDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
drugsDisease_drugsDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DrugsDisease",drugsDisease_drugsDisease_communityDistance);
#drugsChemical-drugsChemical
tmp = g_drugsChemical_drugsChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_drugsChemical_drugsChemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
drugsChemical_drugsChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DrugsChemical",drugsChemical_drugsChemical_communityDistance);
#diseaseChemical-diseaseChemical
tmp = g_diseaseChemical_diseaseChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_diseaseChemical_diseaseChemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
diseaseChemical_diseaseChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_diseaseChemical",diseaseChemical_diseaseChemical_communityDistance);
#nano-drugs-disease
tmp = g_nanoDrugsDisease_nanoDrugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_nanoDrugsDisease_nanoDrugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
nanoDrugsDisease_nanoDrugsDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDrugsDisease",nanoDrugsDisease_nanoDrugsDisease_communityDistance);
#nanoDiseaseChemical
tmp = g_nanoDiseaseChemical_nanoDiseaseChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_nanoDiseaseChemical_nanoDiseaseChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
nanoDiseaseChemical_nanoDiseaseChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_nanoDiseaseChemical",nanoDiseaseChemical_nanoDiseaseChemical_communityDistance);
#drugsDiseaseChemical
tmp = g_drugsDiseaseChemical_drugsDiseaseChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_drugsDiseaseChemical_drugsDiseaseChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
drugsDiseaseChemical_drugsDiseaseChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_drugsDiseaseChemical",drugsDiseaseChemical_drugsDiseaseChemical_communityDistance);
#nanoDrugsChemical
tmp = g_nanoDrugsChemical_nanoDrugsChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_no_all_nano_nanoDrugsChemical_nanoDrugsChemical
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
nanoDrugsChemical_nanoDrugsChemical_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_nanoDrugsChemical",nanoDrugsChemical_nanoDrugsChemical_communityDistance);
#All
tmp = g_W_ADJ_sub
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
tmp = g_W_ADJ_no_all_nano
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
all_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_all",all_communityDistance);
