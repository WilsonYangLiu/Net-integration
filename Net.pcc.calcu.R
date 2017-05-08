# Gave edges with at least two column of tables, this script can calculate the PCC for each edge, based on 
# some measurements (such as, expression data)

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/6.0.expr.cur/cor.Mmp13~Tgfb3/')

#1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
gene.expr <- read.csv(file = '../../0.Backup/[!]fpkm/FPKM_after.txt', sep = '\t', row.names = 1)
sample <- colnames(gene.expr); sample <- grep('X.*', sample, value = TRUE)
gene.expr <- gene.expr[, sample]
#1-----#

# gene.expr <- read.csv(file = '../../2+.Count2FPKM/fpkm.csv', row.names = 1)
# sample <- colnames(gene.expr); sample <- grep('X.*W[^CN].*', sample, value = TRUE)
# gene.expr <- gene.expr[, sample]

save(gene.expr, sample, file = './expr.data.for.cur.RData')

# ---start here--- #
load(file = './expr.data.for.cur.RData')

myfunc <- function(edge.name) {
  Edge <- read.csv(file = paste(edge.name, '.tsv', sep = ''), sep = '\t', stringsAsFactors = FALSE)
  Node <- unique(c(Edge$Node1, Edge$Node2))
  
  tmp <- Node[!(Node %in% stringr::str_to_upper(rownames(gene.expr)))]
  Node.tmp <- gene.expr[stringr::str_to_upper(rownames(gene.expr)) %in% Node, ]
  rownames(Node.tmp) <- stringr::str_to_upper(rownames(Node.tmp))
  
  Edge.tmp <- Edge[!(Edge$Node1 %in% tmp) & !(Edge$Node2 %in% tmp), ]
  # for (ptn in unique(sub('(X.*)WT(.*)', '\\1', x = sample))) {
  #1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
  for (ptn in unique(sub('(X.*)w(.*)', '\\1', x = sample))) {
  #1-----#
    selt <- grep(ptn, sample, value = TRUE)
    # print(selt)
    Edge.tmp[ptn] <- apply(Edge.tmp, MARGIN = 1, FUN = function(x) {
      cor(as.numeric(Node.tmp[x[1], selt]), as.numeric(Node.tmp[x[2], selt]), method = 'pearson')
    })
    # colnames(Edge.tmp) <- c(colnames(Edge.tmp)[1:(dim(Edge.tmp)[2]-1)], ptn)
  }
  x <- Edge[(Edge$Node1 %in% tmp) | (Edge$Node2 %in% tmp), ]
  if (dim(x)[1] > 0) {
    # x$X14 <- 0; x$X17 <- 0; x$X3 <- 0; x$X5 <- 0; x$X9 <- 0
    #1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
    x$X3 <- 0; x$X5 <- 0; x$X9 <- 0; x$X14 <- 0; x$X17 <- 0
    #1-----#
    Edge <- rbind(Edge.tmp, x)
  } else {
    Edge <- Edge.tmp
  }
  Edge[, 4:8] <- round(Edge[, 4:8], digits = 2)
  write.table(Edge, file = paste(edge.name, '.txt', sep = ''), quote = FALSE, sep = '\t', row.names = FALSE)
  
  # for (ptn in unique(sub('(X.*)WT(.*)', '\\1', x = sample))) {
  #1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
  for (ptn in unique(sub('(X.*)w(.*)', '\\1', x = sample))) {
  #1-----#
    selt <- grep(ptn, sample, value = TRUE)
    Node.tmp[ptn] <- rowMeans(Node.tmp[, selt])
    # colnames(Node.tmp) <- c(colnames(Node.tmp)[1:(dim(Node.tmp)[2]-1)], ptn)
  }
  Node.tmp <- Node.tmp[, 26:30]
  x <- as.data.frame(matrix(0, nrow = length(tmp), ncol = 5))
  rownames(x) <- tmp; colnames(x) <- colnames(Node.tmp)
  Node <- rbind(Node.tmp, x)
  Node$Gene <- rownames(Node)
  Node <- Node[, c(6, 1:5)]
  write.table(Node, file = paste(strsplit(edge.name, split = '\\.')[[1]][1], '.node.txt', sep = ''), quote = FALSE, sep = '\t', row.names = FALSE)
}

myfunc('combine.edge')

# Dcn vs Tgfb3
for (ptn in unique(sub('(X.*)WT(.*)', '\\1', x = sample))) {
  selt <- grep(ptn, sample, value = TRUE)
  print(ptn)
  print(cor(as.numeric(gene.expr['Dcn', selt]), as.numeric(gene.expr['Tgfb3', selt]), method = 'pearson'))
}

# Tgfb1 vs Tgfb3
for (ptn in unique(sub('(X.*)w(.*)', '\\1', x = sample))) {
  selt <- grep(ptn, sample, value = TRUE)
  print(ptn)
  print(cor(as.numeric(gene.expr['Tgfb1', selt]), as.numeric(gene.expr['Tgfb3', selt]), method = 'pearson'))
}

# Dcn vs Tgfb1
for (ptn in unique(sub('(X.*)w(.*)', '\\1', x = sample))) {
  selt <- grep(ptn, sample, value = TRUE)
  print(ptn)
  print(cor(as.numeric(gene.expr['Dcn', selt]), as.numeric(gene.expr['Tgfb1', selt]), method = 'pearson'))
}

# Smad2/3 vs Mmp13
for (ptn in unique(sub('(X.*)w(.*)', '\\1', x = sample))) {
  selt <- grep(ptn, sample, value = TRUE)
  print(ptn)
  print(cor(as.numeric(gene.expr['Smad2', selt]), as.numeric(gene.expr['Mmp13', selt]), method = 'pearson'))
}

