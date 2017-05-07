# 

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/tool.network//')

#-- Biogrid
Biogrid <- read.csv(file = './Biogrid/BIOGRID-ORGANISM-Mus_musculus-3.4.144.tab2.txt', sep = '\t')
Biogrid$X.BioGRID.Interaction.ID <- paste('Biogrid', Biogrid$X.BioGRID.Interaction.ID, sep = ':')
Biogrid.flt <- Biogrid[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B", "X.BioGRID.Interaction.ID")]
colnames(Biogrid.flt) <- c('node1', 'node2', 'source'); rm(Biogrid)
Biogrid.flt$node1 <- stringr::str_to_upper(Biogrid.flt$node1)
Biogrid.flt$node2 <- stringr::str_to_upper(Biogrid.flt$node2)
Biogrid.flt$edge <- apply(Biogrid.flt, MARGIN = 1, FUN = function(x) {
  node1 <- x[1]; node2 <- x[2]
  return(paste(sort(c(node1, node2)), collapse = '|'))
})

Biogrid.edge <- unique(Biogrid.flt$edge)
Biogrid.source <- c()
for (item in Biogrid.edge) {
  Biogrid.source <- c(Biogrid.source, paste(Biogrid.flt$source[Biogrid.flt$edge == item], collapse = '|'))
}
Biogrid.flt <- data.frame(node1 = sapply(Biogrid.edge, FUN = function(x) stringr::str_split(x, '\\|')[[1]][1]), 
                          node2 = sapply(Biogrid.edge, FUN = function(x) stringr::str_split(x, '\\|')[[1]][2]),
                          source = Biogrid.source, 
                          stringsAsFactors = FALSE)
save(Biogrid.flt, file = './Biogrid.flt.RData')
write.table(Biogrid.flt, file = gzfile('./Biogrid.flt.tsv.gz', open = 'w'), quote = FALSE, sep = '\t')
#-- [END] Biogrid

#-- String
StringDB <- read.csv(file = gzfile('./StringDB/STRING.edge.mus.tsv.gz'), sep = '\t')
# StringDB$ID <- paste('StringDB', 1:dim(StringDB)[1], sep = ':')
# StringDB <- StringDB[, c('ID', colnames(StringDB)[1:(dim(StringDB)[2]-1)])]
# write.table(StringDB, file = './StringDB/STRING.edge.mus.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
StringDB.flt <- StringDB[(StringDB$experimental != 0) | ((StringDB$experimental == 0) & (StringDB$combined_score >= 700)), 
                         c(2, 3, 1, 4:dim(StringDB)[2])]
colnames(StringDB.flt) <- c('node1', 'node2', 'source'); rm(StringDB)
StringDB.flt$node1 <- stringr::str_to_upper(StringDB.flt$node1)
StringDB.flt$node2 <- stringr::str_to_upper(StringDB.flt$node2)
StringDB.flt$edge <- apply(StringDB.flt, MARGIN = 1, FUN = function(x) {
  node1 <- x[1]; node2 <- x[2]
  return(paste(sort(c(node1, node2)), collapse = '|'))
})

StringDB.edge <- unique(StringDB.flt$edge)
StringDB.source <- c()
for (item in StringDB.edge) {
  StringDB.source <- c(StringDB.source, paste(StringDB.flt$source[StringDB.flt$edge == item], collapse = '|'))
}

tmp <- unique(StringDB.flt[, 4:dim(StringDB.flt)[2]])
rownames(tmp) <- tmp$edge
tmp$node1 <- sapply(tmp$edge, FUN = function(x) stringr::str_split(x, '\\|')[[1]][1])
tmp$node2 <- sapply(tmp$edge, FUN = function(x) stringr::str_split(x, '\\|')[[1]][2])

names(StringDB.source) <- StringDB.edge
tmp$source <- StringDB.source[tmp$edge]
StringDB.flt <- tmp[, c(10:12, 2:9)]; rm(tmp)
save(StringDB.flt, file = './StringDB.flt_expe_0.7.RData')
write.table(StringDB.flt, file = gzfile('./StringDB.flt_expe_0.7.tsv.gz', open = 'w'), quote = FALSE, sep = '\t')
#-- [END] StringDB

# Edges. More stringent for StringDB.flt
StringDB.flt.1 <- StringDB.flt[((StringDB.flt$experimental >= 500) & (StringDB.flt$combined_score >= 700)) | 
                                 ((StringDB.flt$experimental > 0) & (StringDB.flt$experimental < 500) & (StringDB.flt$combined_score >= 900)) | 
                                 ((StringDB.flt$experimental == 0) & (StringDB.flt$combined_score >= 950)), 1:3]
tmp <- intersect(rownames(StringDB.flt.1), rownames(Biogrid.flt))
StringDB.flt.1[tmp, 3] <- paste(StringDB.flt.1[tmp, 3], Biogrid.flt[tmp, 3], sep = '|')
Edges <- rbind(StringDB.flt.1, Biogrid.flt[!(rownames(Biogrid.flt) %in% tmp), ])
write.table(Edges, file = './Edge.tsv', quote = FALSE, sep = '\t')

# Nodes
Nodes <- unique(c(Edges$node1, Edges$node2))
Nodes <- data.frame(node = Nodes, 
                    row.names = Nodes, 
                    stringsAsFactors = FALSE)
Nodes$isRegulator <- '-'
Nodes$AnimalTFDB <- '-'

#-- AnimalTFDB
load(file = './AnimalTFDB/Mus_musculus.RData')
for(n in names(MUS)) {
  MUS[[n]] <- sapply(MUS[[n]], FUN = function(x) stringr::str_to_upper(x))
}
Nodes[Nodes$node %in% MUS$TF, 'AnimalTFDB'] <- 'TF'
Nodes[Nodes$node %in% MUS$TcF, 'AnimalTFDB'] <- 'TcF'
Nodes[Nodes$node %in% MUS$CRF, 'AnimalTFDB'] <- 'CRF'
#-- [END] AnimalTFDB

Nodes$RegNetwork <- '-'
#-- RegNetwork
RegNetwork <- read.csv(file = './RegNetwork/mouse.source', sep = '\t', stringsAsFactors = FALSE)
RegNode <- sapply(RegNetwork$Regulator.Symbol, FUN = function(x) stringr::str_to_upper(x))
RegNode <- unique(RegNode)
Nodes[Nodes$node %in% RegNode, 'RegNetwork'] <- '+'
#-- [END] RegNetwork

# Nodes
Nodes[apply(Nodes, MARGIN = 1, FUN = function(x) {(x[3] != '-') | (x[4] != '-')}), 'isRegulator'] <- '+'
write.table(Nodes, file = './Node.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

