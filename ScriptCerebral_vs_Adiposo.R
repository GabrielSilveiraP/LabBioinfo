#Primeiro script tentando buscar a expressao diferencial
library(DESeq2)
library(tidyverse)

tabelagenes <- read.csv("C:/UFSC/LabInfo/prepDE/tabelagenes.csv", row.names = 1 )
head(tabelagenes)
tabelatranscritos <- read.csv("C:/UFSC/LabInfo/prepDE/tabelatranscritos.csv", row.names =  1 )
head(tabelatranscritos)
#Porq os dois colnames e n sendo um row? Por causa que os sao somente colunas, vai ter diferenças nas linhas
all(colnames(tabelagenes) %in% colnames(tabelatranscritos))
#Output -> true, ou seja, tem as mesmas amostras
all(colnames(tabelagenes) == colnames(tabelatranscritos))
#O de cima é p ver se eles estao na mesma ordem

DadosAmostras <- data.frame(
    amostras = c("SRR21831518", "SRR21831519", "SRR21831520", 
    "SRR21837144", "SRR21837145", "SRR21837146"),
    tecido = c("cerebral", "cerebral", "cerebral",
                "adiposo", "adiposo", "adiposo")
)
  
rownames(DadosAmostras) <-  DadosAmostras$amostras
#P transformar o frame ali em um pedaço da tabela
all(colnames(tabelagenes) == rownames(DadosAmostras))
#Verificar a ordem das amostras entre a tabgene e amostras

#Juntando os dados com o DEseq2
dds <-DESeqDataSetFromMatrix(
  countData = tabelatranscritos,
  colData = DadosAmostras,
  design =  ~ tecido)

dds

#Filtragem p retirar alguns reads muito pequenos
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep]

dds

#setando A referencia(Pq cerebral? Porq o rangeli n tem historico de infectar cerebral, \
#ent estou botando aqlq q eu espero que tenha menos expressao )
dds$tecido <- relevel(dds$tecido, ref = "adiposo")
#Se n fizer isso aq, ele vai usa ro primeiro em ordem alfabetica

#Rodando o deseqdds
dds <- DESeq(dds)

res <- results(dds)

res


#ja tamos nos resultados aqui

summary(res)
#Fazendo com q o p value seja de 0.01 e n 0.1
res0.01 <-  results(dds,alpha = 0.01)

summary(res0.01)


#Plottando

plotMA(res)
plotMA(res0.01)
