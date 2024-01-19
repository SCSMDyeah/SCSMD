select_method <-function(labels_SC3,labels_Seurat,labels_SHARP,labels_cidr,labels_SINCERA,labels_Rphenograph,labels_RaceID){



#sc3
sc31 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_Seurat))
sc32 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_SHARP))
sc33 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_cidr))
sc34 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_SINCERA))
sc35 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_Rphenograph))
sc36 <- adjustedRandIndex(as.vector(labels_SC3),as.vector(labels_RaceID))
sc3_AVG <- sum(c(sc31,sc32,sc33,sc34,sc35,sc36))/6


#Seurat
seurtar1 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_SC3))
seurtar2 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_SHARP))
seurtar3 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_cidr))
seurtar4 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_SINCERA))
seurtar5 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_Rphenograph))
seurtar6 <- adjustedRandIndex(as.vector(labels_Seurat),as.vector(labels_RaceID))
seurtar_AVG <- sum(c(seurtar1,seurtar2,seurtar3,seurtar4,seurtar5,seurtar6))/6


#SHARP
SHARP1 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_SC3))
SHARP2 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_Seurat))
SHARP3 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_cidr))
SHARP4 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_SINCERA))
SHARP5 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_Rphenograph))
SHARP6 <- adjustedRandIndex(as.vector(labels_SHARP),as.vector(labels_RaceID))
SHARP_AVG <- sum(c(SHARP1,SHARP2,SHARP3,SHARP4,SHARP5,SHARP6))/6


#cidr
cidr1 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_SC3))
cidr2 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_Seurat))
cidr3 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_SHARP))
cidr4 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_SINCERA))
cidr5 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_Rphenograph))
cidr6 <- adjustedRandIndex(as.vector(labels_cidr),as.vector(labels_RaceID))
cidr_AVG <- sum(c(cidr1,cidr2,cidr3,cidr4,cidr5,cidr6))/6


#SINCERA
SINCERA1 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_SC3))
SINCERA2 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_Seurat))
SINCERA3 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_SHARP))
SINCERA4 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_cidr))
SINCERA5 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_Rphenograph))
SINCERA6 <- adjustedRandIndex(as.vector(labels_SINCERA),as.vector(labels_RaceID))
SINCERA_AVG <- sum(c(SINCERA1,SINCERA2,SINCERA3,SINCERA4,SINCERA5,SINCERA6))/6


#Rphenograph
Rphenograph1 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_SC3))
Rphenograph2 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_Seurat))
Rphenograph3 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_SHARP))
Rphenograph4 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_cidr))
Rphenograph5 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_SINCERA))
Rphenograph6 <- adjustedRandIndex(as.vector(labels_Rphenograph),as.vector(labels_RaceID))
Rphenograph_AVG <- sum(c(Rphenograph1,Rphenograph2,Rphenograph3,Rphenograph4,Rphenograph5,Rphenograph6))/6


#RaceID
RaceID1 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_SC3))
RaceID2 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_Seurat))
RaceID3 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_SHARP))
RaceID4 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_cidr))
RaceID5 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_SINCERA))
RaceID6 <- adjustedRandIndex(as.vector(labels_RaceID),as.vector(labels_Rphenograph))
RaceID_AVG <- sum(c(RaceID1,RaceID2,RaceID3,RaceID4,RaceID5,RaceID6))/6

ran <-matrix(c(sc3_AVG, seurtar_AVG, SHARP_AVG, cidr_AVG, SINCERA_AVG, Rphenograph_AVG, RaceID_AVG)) 
rownames(ran) <-c("sc3_AVG", "seurat_AVG", "SHARP_AVG", "cidr_AVG", "SINCERA_AVG", "Rphenograph_AVG", "RaceID_AVG")
method_rank <-matrix(rank(ran))   #4 2 6 3 7 1 5
rownames(method_rank) <-c("sc3", "seurat", "SHARP", "cidr", "SINCERA", "Rphenograph", "RaceID")

select <- list(ran, method_rank)
names(select)[1] <- "method_AVG"
names(select)[2] <- "rank"
return(select)
}