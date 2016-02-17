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
load("~/unisa/II_anno/Bioinformatica/progetto/myWork/.RData")

#ADJ = distance matrix without CLR 
#SubMatrix without chemical row-col
ADJ_no_chemical = rbind(cbind(ADJ[nano,nano],
                        ADJ[nano,drugs],
                        ADJ[nano,disease]),
                        cbind(ADJ[drugs,nano],
                              ADJ[drugs,drugs],
                              ADJ[drugs,disease]),
                        cbind(ADJ[disease,nano],
                              ADJ[disease,drugs],
                              ADJ[disease,disease]));
#dump("KTI", file='KTI.r');
#Computing CLR
W_ADJ = KTI(KTDM_mat = abs(ADJ),SIGN = sign(ADJ))
W_ADJ_sub = rbind(cbind(W_ADJ[nano,nano],
                        W_ADJ[nano,drugs],
                        W_ADJ[nano,disease]),
                  cbind(W_ADJ[drugs,nano],
                        W_ADJ[drugs,drugs],
                        W_ADJ[drugs,disease]),
                  cbind(W_ADJ[disease,nano],
                        W_ADJ[disease,drugs],
                        W_ADJ[disease,disease]));
W_ADJ_no_chemical = KTI(KTDM_mat = abs(ADJ_no_chemical),SIGN = sign(ADJ_no_chemical))

#Getting submatrices to compare
W_ADJ_nano_nano = W_ADJ[nano,nano];
W_ADJ_nano_drugs = W_ADJ[nano,drugs];
W_ADJ_nano_disease = W_ADJ[nano,disease];
W_ADJ_drugs_nano = W_ADJ[drugs,nano];
W_ADJ_drugs_drugs = W_ADJ[drugs,drugs];
W_ADJ_drugs_disease = W_ADJ[drugs,disease];
W_ADJ_disease_nano = W_ADJ[disease, nano];
W_ADJ_disease_drugs = W_ADJ[disease, drugs];
W_ADJ_disease_disease = W_ADJ[disease, disease];
W_ADJ_nanoDrugs_nanoDrugs = rbind(
  cbind(W_ADJ[nano,nano],W_ADJ[nano,drugs]),
  cbind(W_ADJ[drugs,nano],W_ADJ[drugs,drugs])
);
W_ADJ_nanoDisease_nanoDisease = rbind(
  cbind(W_ADJ[nano,nano],W_ADJ[nano,disease]),
  cbind(W_ADJ[disease,nano],W_ADJ[disease,disease])
);
W_ADJ_drugsDisease_drugsDisease = rbind(
  cbind(W_ADJ[drugs,drugs],W_ADJ[drugs,disease]),
  cbind(W_ADJ[disease,drugs],W_ADJ[disease,disease])
);

W_ADJ_no_chemical_nano_nano = W_ADJ_no_chemical[nano,nano];
W_ADJ_no_chemical_nano_drugs = W_ADJ_no_chemical[nano,drugs];
W_ADJ_no_chemical_nano_disease = W_ADJ_no_chemical[nano,disease];
W_ADJ_no_chemical_drugs_nano = W_ADJ_no_chemical[drugs,nano];
W_ADJ_no_chemical_drugs_drugs = W_ADJ_no_chemical[drugs,drugs];
W_ADJ_no_chemical_drugs_disease = W_ADJ_no_chemical[drugs,disease];
W_ADJ_no_chemical_disease_nano = W_ADJ_no_chemical[disease, nano];
W_ADJ_no_chemical_disease_drugs = W_ADJ_no_chemical[disease, drugs];
W_ADJ_no_chemical_disease_disease = W_ADJ_no_chemical[disease, disease];
W_ADJ_no_chemical_nanoDrugs_nanoDrugs = rbind(
  cbind(W_ADJ_no_chemical[nano,nano],W_ADJ_no_chemical[nano,drugs]),
  cbind(W_ADJ_no_chemical[drugs,nano],W_ADJ_no_chemical[drugs,drugs])
);
W_ADJ_no_chemical_nanoDisease_nanoDisease = rbind(
  cbind(W_ADJ_no_chemical[nano,nano],W_ADJ_no_chemical[nano,disease]),
  cbind(W_ADJ_no_chemical[disease,nano],W_ADJ_no_chemical[disease,disease])
);
W_ADJ_no_chemical_drugsDisease_drugsDisease = rbind(
  cbind(W_ADJ_no_chemical[drugs,drugs],W_ADJ_no_chemical[drugs,disease]),
  cbind(W_ADJ_no_chemical[disease,drugs],W_ADJ_no_chemical[disease,disease])
);

#Nano-Nano analysis plot
tmp = W_ADJ_nano_nano;
diag(tmp) = 0;
g = graph.adjacency(tmp,mode="undirected",weighted = T)
plot(g,edge.width = E(g)$weight * 1.2)

dev.new()
tmp = W_ADJ_no_chemical_nano_nano;
diag(tmp) = 0;
g1 = graph.adjacency(tmp,mode="undirected",weighted = T)
plot(g1,edge.width = E(g1)$weight * 1.2)

par(mfrow=c(1,2))
plot(g,edge.width = E(g)$weight * 1.2, main="Nano_nano complete CLR")
plot(g1,edge.width = E(g1)$weight * 1.2, main="Nano_nano CLR without chemical")

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
cluster_W_ADJ_no_chemical_nano_nano <- kmeans(W_ADJ_no_chemical_nano_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_nano$cluster,cluster_W_ADJ_no_chemical_nano_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot <- ggplot(confusion)
plot + ggtitle("Nano-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-drugs
cluster_W_ADJ_nano_drugs <- kmeans(W_ADJ_nano_drugs,k,iterations)
cluster_W_ADJ_no_chemical_nano_drugs <- kmeans(W_ADJ_no_chemical_nano_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_drugs$cluster,cluster_W_ADJ_no_chemical_nano_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Nano-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-disease
cluster_W_ADJ_nano_disease <- kmeans(W_ADJ_nano_disease,k,iterations)
cluster_W_ADJ_no_chemical_nano_disease <- kmeans(W_ADJ_no_chemical_nano_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nano_disease$cluster,cluster_W_ADJ_no_chemical_nano_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Nano-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-nano
cluster_W_ADJ_drugs_nano <- kmeans(W_ADJ_drugs_nano,k,iterations)
cluster_W_ADJ_no_chemical_drugs_nano <- kmeans(W_ADJ_no_chemical_drugs_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_nano$cluster,cluster_W_ADJ_no_chemical_drugs_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-drugs
cluster_W_ADJ_drugs_drugs <- kmeans(W_ADJ_drugs_drugs,k,iterations)
cluster_W_ADJ_no_chemical_drugs_drugs <- kmeans(W_ADJ_no_chemical_drugs_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_drugs$cluster,cluster_W_ADJ_no_chemical_drugs_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugs-disease
cluster_W_ADJ_drugs_disease <- kmeans(W_ADJ_drugs_disease,k,iterations)
cluster_W_ADJ_no_chemical_drugs_disease <- kmeans(W_ADJ_no_chemical_drugs_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugs_disease$cluster,cluster_W_ADJ_no_chemical_drugs_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Drugs-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-nano
cluster_W_ADJ_disease_nano <- kmeans(W_ADJ_disease_nano,k,iterations)
cluster_W_ADJ_no_chemical_disease_nano <- kmeans(W_ADJ_no_chemical_disease_nano,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_nano$cluster,cluster_W_ADJ_no_chemical_disease_nano$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Nano clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-drugs
cluster_W_ADJ_disease_drugs <- kmeans(W_ADJ_disease_drugs,k,iterations)
cluster_W_ADJ_no_chemical_disease_drugs <- kmeans(W_ADJ_no_chemical_disease_drugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_drugs$cluster,cluster_W_ADJ_no_chemical_disease_drugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Drugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#disease-disease
cluster_W_ADJ_disease_disease <- kmeans(W_ADJ_disease_disease,k,iterations)
cluster_W_ADJ_no_chemical_disease_disease <- kmeans(W_ADJ_no_chemical_disease_disease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_disease_disease$cluster,cluster_W_ADJ_no_chemical_disease_disease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("Disease-Disease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nanoDrugs-nanoDrugs
cluster_W_ADJ_nanoDrugs_nanoDrugs <- kmeans(W_ADJ_nanoDrugs_nanoDrugs,k,iterations)
cluster_W_ADJ_no_chemical_nanoDrugs_nanoDrugs <- kmeans(W_ADJ_no_chemical_nanoDrugs_nanoDrugs,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDrugs_nanoDrugs$cluster,cluster_W_ADJ_no_chemical_nanoDrugs_nanoDrugs$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDrugs-NanoDrugs clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nanoDisease-nanoDisease
cluster_W_ADJ_nanoDisease_nanoDisease <- kmeans(W_ADJ_nanoDisease_nanoDisease,k,iterations)
cluster_W_ADJ_no_chemical_nanoDisease_nanoDisease <- kmeans(W_ADJ_no_chemical_nanoDisease_nanoDisease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_nanoDisease_nanoDisease$cluster,cluster_W_ADJ_no_chemical_nanoDisease_nanoDisease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDisease-NanoDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#drugsDisease-drugsDisease
cluster_W_ADJ_drugsDisease_drugsDisease <- kmeans(W_ADJ_drugsDisease_drugsDisease,k,iterations)
cluster_W_ADJ_no_chemical_drugsDisease_drugsDisease <- kmeans(W_ADJ_no_chemical_drugsDisease_drugsDisease,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ_drugsDisease_drugsDisease$cluster,cluster_W_ADJ_no_chemical_drugsDisease_drugsDisease$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("DrugsDisease-DrugsDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 
#nano-drugs-disease
cluster_W_ADJ <- kmeans(W_ADJ_sub,k,iterations)
cluster_W_ADJ_no_chemical <- kmeans(W_ADJ_no_chemical,k,iterations);
tmp = confusionMatrix(table(cluster_W_ADJ$cluster,cluster_W_ADJ_no_chemical$cluster))
confusion <- as.data.frame(as.table(tmp$table));
plot #restore plot
plot <- ggplot(confusion)
plot + ggtitle("NanoDrugsDisease-NanoDrugsDisease clustering") + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="I result") + scale_y_discrete(name="II result") 

#Social network analysis test
#Nano-Nano analysis
g_nano_nano = graph.adjacency(W_ADJ_nano_nano,mode="undirected",weighted = T)
g_no_chemical_nano_nano = graph.adjacency(W_ADJ_no_chemical_nano_nano,mode="undirected",weighted = T)
degree_nano_nano = degree(g_nano_nano);
degree_no_chemical_nano_nano = degree(g_no_chemical_nano_nano);
tmp=sort(degree_nano_nano,decreasing = TRUE);
tmp1=sort(degree_no_chemical_nano_nano,decreasing = TRUE);
degreeNanoNano_distance = kendallTauDistance(names(tmp), names(tmp1));
table = rbind(cbind('KendallTau Distance',degreeNanoNano_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,'degree_nanoNano.pdf');
# plot(degree_nano_nano, degree_no_chemical_nano_nano);
# plot(g_nano_nano,
#      vertex.size = round(degree_nano_nano, 2),
#      vertex.color = rainbow(vcount(g_nano_nano)))
# dev.new()
# plot(g_no_chemical_nano_nano,
#      vertex.size = round(degree_no_chemical_nano_nano, 2),
#      vertex.color = rainbow(vcount(g_no_chemical_nano_nano)))
#Drugs-Drugs analysis
g_drugs_drugs = graph.adjacency(W_ADJ_drugs_drugs,mode="undirected",weighted = T)
g_no_chemical_drugs_drugs = graph.adjacency(W_ADJ_no_chemical_drugs_drugs,mode="undirected",weighted = T)
degree_drugs_drugs = degree(g_drugs_drugs);
degree_no_chemical_drugs_drugs = degree(g_no_chemical_drugs_drugs);
tmp=sort(degree_drugs_drugs,decreasing = TRUE);
tmp1=sort(degree_no_chemical_drugs_drugs,decreasing = TRUE);
degreeDrugsDrugs_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsDrugs_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_drugsDrugs.pdf");
#Disease-Disease analysis
g_disease_disease = graph.adjacency(W_ADJ_disease_disease,mode="undirected",weighted = T)
g_no_chemical_disease_disease = graph.adjacency(W_ADJ_no_chemical_disease_disease,mode="undirected",weighted = T)
degree_disease_disease = degree(g_disease_disease);
degree_no_chemical_disease_disease = degree(g_no_chemical_disease_disease);
tmp=sort(degree_disease_disease,decreasing = TRUE);
tmp1=sort(degree_no_chemical_disease_disease,decreasing = TRUE);
degreeDiseaseDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDiseaseDisease_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_diseaseDisease.pdf");
#NanoDrugs-NanoDrugs analysis
g_nanoDrugs_nanoDrugs = graph.adjacency(W_ADJ_nanoDrugs_nanoDrugs,mode="undirected",weighted = T)
g_no_chemical_nanoDrugs_nanoDrugs = graph.adjacency(W_ADJ_no_chemical_nanoDrugs_nanoDrugs,mode="undirected",weighted = T)
degree_nanoDrugs_nanoDrugs = degree(g_nanoDrugs_nanoDrugs);
degree_no_chemical_nanoDrugs_nanoDrugs = degree(g_no_chemical_nanoDrugs_nanoDrugs);
tmp=sort(degree_nanoDrugs_nanoDrugs,decreasing = TRUE);
tmp1=sort(degree_no_chemical_nanoDrugs_nanoDrugs,decreasing = TRUE);
degreeNanoDrugs_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDrugs_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_nanoDrugs_nanoDrugs.pdf");
#NanoDisease-NanoDisease analysis
g_nanoDisease_nanoDisease = graph.adjacency(W_ADJ_nanoDisease_nanoDisease,mode="undirected",weighted = T)
g_no_chemical_nanoDisease_nanoDisease = graph.adjacency(W_ADJ_no_chemical_nanoDisease_nanoDisease,mode="undirected",weighted = T)
degree_nanoDisease_nanoDisease = degree(g_nanoDisease_nanoDisease);
degree_no_chemical_nanoDisease_nanoDisease = degree(g_no_chemical_nanoDisease_nanoDisease);
tmp=sort(degree_nanoDisease_nanoDisease,decreasing = TRUE);
tmp1=sort(degree_no_chemical_nanoDisease_nanoDisease,decreasing = TRUE);
degreeNanoDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeNanoDisease_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_nanoDisease_nanoDisease.pdf");
#DrugDisease-DrugDisease analysis
g_drugsDisease_drugsDisease = graph.adjacency(W_ADJ_drugsDisease_drugsDisease,mode="undirected",weighted = T)
g_no_chemical_drugsDisease_drugsDisease = graph.adjacency(W_ADJ_no_chemical_drugsDisease_drugsDisease,mode="undirected",weighted = T)
degree_drugsDisease_drugsDisease = degree(g_drugsDisease_drugsDisease);
degree_no_chemical_drugsDisease_drugsDisease = degree(g_no_chemical_drugsDisease_drugsDisease);
tmp=sort(degree_drugsDisease_drugsDisease,decreasing = TRUE);
tmp1=sort(degree_no_chemical_drugsDisease_drugsDisease,decreasing = TRUE);
degreeDrugsDisease_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeDrugsDisease_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_drugsDisease_drugsDisease.pdf");
#all matrix
g_W_ADJ_no_chemical = graph.adjacency(W_ADJ_no_chemical,mode="undirected",weighted = T)
g_W_ADJ_sub = graph.adjacency(W_ADJ_sub,mode="undirected",weighted = T)
degree_W_ADJ_sub = degree(g_W_ADJ_sub);
degree_W_ADJ_no_chemical = degree(g_W_ADJ_no_chemical);
tmp=sort(degree_W_ADJ_sub,decreasing = TRUE);
tmp1=sort(degree_W_ADJ_no_chemical,decreasing = TRUE);
degreeAll_distance = kendallTauDistance(tmp,tmp1);
table = rbind(cbind('KendallTau Distance',degreeAll_distance,'',''),cbind('DEGREE (no chemical)','NAME','DEGREE (with chemical)','NAME'),cbind(as.numeric(tmp),names(tmp),as.numeric(tmp1),names(tmp1)));
plotTable(table,"degree_All.pdf");

#Communities check
#nano-nano
par(mfrow=c(1,2))
plot(walktrap.community(g_nano_nano),simplify(g_nano_nano,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_nano_nano),simplify(g_no_chemical_nano_nano,remove.loops = TRUE))
#drugs-drugs
par(mfrow=c(1,2))
plot(walktrap.community(g_drugs_drugs),simplify(g_drugs_drugs,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_drugs_drugs),simplify(g_no_chemical_drugs_drugs,remove.loops = TRUE))
#disease-disease
par(mfrow=c(1,2))
plot(walktrap.community(g_disease_disease),simplify(g_disease_disease,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_disease_disease),simplify(g_no_chemical_disease_disease,remove.loops = TRUE))
#nanoDrugs-nanoDrugs
par(mfrow=c(1,2))
plot(walktrap.community(g_nanoDrugs_nanoDrugs),simplify(g_nanoDrugs_nanoDrugs,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_nanoDrugs_nanoDrugs),simplify(g_no_chemical_nanoDrugs_nanoDrugs,remove.loops = TRUE))
#nanoDisease-nanoDisease
par(mfrow=c(1,2))
plot(walktrap.community(g_nanoDisease_nanoDisease),simplify(g_nanoDisease_nanoDisease,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_nanoDisease_nanoDisease),simplify(g_no_chemical_nanoDisease_nanoDisease,remove.loops = TRUE))
#drugsDisease-drugsDisease
par(mfrow=c(1,2))
plot(walktrap.community(g_drugsDisease_drugsDisease),simplify(g_drugsDisease_drugsDisease,remove.loops = TRUE))
plot(walktrap.community(g_no_chemical_drugsDisease_drugsDisease),simplify(g_no_chemical_drugsDisease_drugsDisease,remove.loops = TRUE))
#nano-drugs-disease
par(mfrow=c(1,2))
plot(walktrap.community(g_W_ADJ_sub),simplify(g_W_ADJ_sub,remove.loops = TRUE))
plot(walktrap.community(g_W_ADJ_no_chemical),simplify(g_W_ADJ_no_chemical,remove.loops = TRUE))

#Interactive plot
#nano-nano
tmp = g_nano_nano;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
edgeColor = matrix(0,dim(W_ADJ_nano_nano)[1],dim(W_ADJ_nano_nano)[2]);
edgeColor[W_ADJ_nano_nano<0]='red';
edgeColor[W_ADJ_nano_nano>=0]='green';
tkplot(simplify(g_nano_nano,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)], edge.color=edgeColor,edge.curved =TRUE)
#gd <- get.data.frame(g_nano_nano, what = "edges")
#simpleNetwork(gd, fontSize = 12)
tmp = g_no_chemical_nano_nano;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
edgeColor = matrix(0,dim(W_ADJ_no_chemical_nano_nano)[1],dim(W_ADJ_no_chemical_nano_nano)[2]);
edgeColor[W_ADJ_no_chemical_nano_nano<0]='red';
edgeColor[W_ADJ_no_chemical_nano_nano>=0]='green';
tkplot(simplify(g_no_chemical_nano_nano,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)], edge.color=edgeColor)
nano_nano_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoNano",nano_nano_communityDistance);
#drugs-drugs
tmp = g_drugs_drugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_drugs_drugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_chemical_drugs_drugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_chemical_drugs_drugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
drugs_drugs_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DrugsDrugs",drugs_drugs_communityDistance);
#disease-disease
tmp = g_disease_disease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_disease_disease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_chemical_disease_disease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_chemical_disease_disease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
disease_disease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DiseaseDisease",disease_disease_communityDistance);
#nanoDrugs-nanoDrugs
tmp = g_nanoDrugs_nanoDrugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_nanoDrugs_nanoDrugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_chemical_nanoDrugs_nanoDrugs;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_chemical_nanoDrugs_nanoDrugs,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
nanoDrugs_nanoDrugs_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDrugs",nanoDrugs_nanoDrugs_communityDistance);
#nanoDisease-nanoDisease
tmp = g_nanoDisease_nanoDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_nanoDisease_nanoDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_chemical_nanoDisease_nanoDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_chemical_nanoDisease_nanoDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
nanoDisease_nanoDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDisease",nanoDisease_nanoDisease_communityDistance);
#drugsDisease-drugsDisease
tmp = g_drugsDisease_drugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_drugsDisease_drugsDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_no_chemical_drugsDisease_drugsDisease;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_no_chemical_drugsDisease_drugsDisease,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
drugsDisease_drugsDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_DrugsDisease",drugsDisease_drugsDisease_communityDistance);
#nano-drugs-disease
tmp = g_W_ADJ_sub;
E(tmp)$weight = abs(E(tmp)$weight)
comm1 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm1))
#tkplot(simplify(g_W_ADJ_sub,remove.loops = TRUE), vertex.color=mycolors[membership(comm1)])
tmp = g_W_ADJ_no_chemical;
E(tmp)$weight = abs(E(tmp)$weight)
comm2 <- fastgreedy.community(tmp);
mycolors <- heat.colors(length(comm2))
#tkplot(simplify(g_W_ADJ_no_chemical,remove.loops = TRUE), vertex.color=mycolors[membership(comm2)])
nanoDrugsDisease_nanoDrugsDisease_communityDistance = compare(comm1, comm2,method='vi')
plotCommunities(comm1,comm2,"Communities_NanoDrugsDisease",nanoDrugsDisease_nanoDrugsDisease_communityDistance);

