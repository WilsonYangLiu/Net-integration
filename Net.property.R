# For a network without any properties in it, we need to find out what properties it contains. 
# This is how this script works. Say we have measurement (such as expression values) for every 
# nodes in the network, we calculate the SD for nodes and the PCC for edges, compared with control.
# NOTE: The `maps` is how the controls should be chosen for a sepecific treatment, say, 
#       w3 to w3_c, w5 to w5_c (w* is the treatment sample, w*_c is the coresponding control sample)

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/')

# Load Global Network
Edges <- read.csv(file = './5.0.network/Edge.mus.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
Nodes <- read.csv(file = './5.0.network/Node.mus.tsv', sep = '\t', stringsAsFactors = FALSE)

# Load Expresion Value
#1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
Gene.count <- read.csv(file = './0.Backup/[!]fpkm/FPKM_after.txt', sep = '\t', row.names = 1)
rownames(Gene.count) <- sapply(rownames(Gene.count), FUN = function(x) stringr::str_to_upper(x))
#1-----#

# Gene.count <- read.csv(file = './2+.Count2FPKM/fpkm.csv', row.names = 1)
# Gene.count <- read.csv(file = './2.count/htseq/count.txt', sep = '\t', row.names = 1)
# rownames(Gene.count) <- sapply(rownames(Gene.count), FUN = function(x) stringr::str_to_upper(x))

# Filter Network
Edges <- Edges[apply(Edges[, c('node1', 'node2')], MARGIN = 1, FUN = function(x) {
  sum(x %in% rownames(Gene.count)) == 2
}), ]
Nodes <- Nodes[Nodes$node %in% unique(c(Edges$node1, Edges$node2)), ]

# All Genes, focus on DEG and DNB
DNB <- read.table(file = './3.2.dnb/module_gene.txt', stringsAsFactors = FALSE)$V1
DNB <- sapply(DNB, FUN = function(x) stringr::str_to_upper(x))
DNB.in.net <- names(DNB[DNB %in% Nodes$node])
write.table(DNB.in.net, file = './5.0.network/DNB.in.flt.network.txt', quote = FALSE, sep = '\n', col.names = FALSE, row.names = FALSE)
Nodes$isDNB <- '-'
Nodes$isDNB[Nodes$node %in% DNB] <- '+'

# Sperate Expression value into differ time points
#1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
wx.count <- list() # For disease 
wx.count.c <- list() # For control
wx.count$w3 <- Gene.count[Nodes$node, grep(pattern = 'X3', x = colnames(Gene.count), value = TRUE)]
wx.count$w5 <- Gene.count[Nodes$node, grep(pattern = 'X5', x = colnames(Gene.count), value = TRUE)]
wx.count$w9 <- Gene.count[Nodes$node, grep(pattern = 'X9', x = colnames(Gene.count), value = TRUE)]
wx.count$w14 <- Gene.count[Nodes$node, grep(pattern = 'X14', x = colnames(Gene.count), value = TRUE)]
wx.count$w17 <- Gene.count[Nodes$node, grep(pattern = 'X17', x = colnames(Gene.count), value = TRUE)]
wx.count.c$wx <- Gene.count[Nodes$node, grep(pattern = 'c', x = colnames(Gene.count), value = TRUE)]
#1-----#

# wx.count <- list() # For disease
# wx.count.c <- list() # For control
# wx.count$w3 <- Gene.count[Nodes$node, grep(pattern = 'X3WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w5 <- Gene.count[Nodes$node, grep(pattern = 'X5WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w9 <- Gene.count[Nodes$node, grep(pattern = 'X9WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w14 <- Gene.count[Nodes$node, grep(pattern = 'X14WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w17 <- Gene.count[Nodes$node, grep(pattern = 'X17WT', x = colnames(Gene.count), value = TRUE)]
# wx.count.c$wx <- Gene.count[Nodes$node, grep(pattern = 'N|C', x = colnames(Gene.count), value = TRUE)]

# mapping between disease and control
maps <- rep('wx', 5); names(maps) <- names(wx.count)

#-- for week x (3, 5, 9, 14, 17)
for (w in names(wx.count)) {
  wx <- wx.count[[w]]
  Nodes[w] <- sapply(Nodes$node, FUN = function(x) {
    sd(wx[x, ])
  })
  
  wx <- t(wx)
  tmp <- apply(Edges[, c('node1', 'node2')], MARGIN = 1, FUN = function(x) {
    cor(x = wx[, x[1]], y = wx[, x[2]], method = "spearman")
  })
  Edges[w] <- tmp
}

for (w in names(wx.count.c)) {
  wx <- wx.count.c[[w]]
  Nodes[w] <- sapply(Nodes$node, FUN = function(x) {
    sd(wx[x, ])
  })
  
  wx <- t(wx)
  tmp <- apply(Edges[, c('node1', 'node2')], MARGIN = 1, FUN = function(x) {
    cor(x = wx[, x[1]], y = wx[, x[2]], method = "spearman")
  })
  Edges[w] <- tmp
}

for (m in 1:length(maps)) {
  Nodes[paste(names(maps[m]), 'adj', sep = '_')] <- Nodes[names(maps[m])] / (Nodes[maps[m]] + 0.001)
  
  Edges[paste(names(maps[m]), 'adj', sep = '_')] <- abs(Edges[names(maps[m])]) / (abs(Edges[maps[m]]) + 0.001)
}

# for (w in names(wx.count.c)) {
#   Nodes[w] <- NULL
#   Edges[w] <- NULL
# }

# Save the Network
write.csv(Nodes, file = './5.0.network/Node.mus.flt.csv', row.names = FALSE, quote = FALSE)
write.csv(Edges, file = './5.0.network/Edge.mus.flt.csv', quote = FALSE)

# Determine the range
max(apply(Nodes[, grep('adj', colnames(Nodes), value = TRUE)], MARGIN = 2, FUN = function(x) max(x)))
tmp <- c(Nodes$w3_adj, Nodes$w5_adj, Nodes$w9_adj, Nodes$w14_adj, Nodes$w17_adj) # 0, 2, 570.0757

max(apply(Edges[, grep('adj', colnames(Nodes), value = TRUE)], MARGIN = 2, FUN = function(x) max(x)))
tmp <- c(Edges$w3_adj, Edges$w5_adj, Edges$w9_adj, Edges$w14_adj, Edges$w17_adj) # 0, 10, 1000
