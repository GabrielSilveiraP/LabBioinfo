library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(tximport)
library(txdbmaker)
library(GenomeInfoDbData)
#-----1 -> Só jogando as coisas p trabaiar
diretorio_base <-"/home/gabriel/Documents/CoisasdoLab"
#Aqui já vai ser diferente porq os arquivos tao escondidos dentro de outras pastas
arquivos_amostras <- list.files(diretorio_base, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
#list files é um comando para que ele ache cada arqv com determinados *patterns*
names(arquivos_amostras) <- basename(dirname(arquivos_amostras))
#O de cima da o nome da pasta para o caminho, é meio confuso mesmo
#tbm, o dirname(DIRETORIO nome) pega o caminho da pasta, basename pega só o nome dela que geralmente vai ser o
#id dela
arquivos_amostras
#checando
#--------2 -> Começando a formar as análises
#aqui q a parada de transcripto vira gene 
gtf_file <- "/home/gabriel/Documents/CoisasdoLab/Trangeli_SC58_84v2.gtf.gz"
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf", organism = "Trypanosoma rangeli")
#importando o de ref
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
head(tx2gene)
#Ai ficou o esquema certinho p puxar as coisas, meio q criamos o molde e agr vamos encher com as inf
#----------2,5 arrumando que os ids n tao batendo
quant_teste <- read_tsv(arquivos_amostras[1])
head(quant_teste$Name)

gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)
colnames(gtf_df)
head(gtf_df)
#Foi p ver como que tava cada nome e achar a diferença
todosIDS <- quant_teste$Name
tx2gene_salmon <- data.frame(
  TXNAME = todosIDS,
  GENEID = sub("\\.t\\d+$", "", todosIDS)  # remove .t1, .t2, que é o final da id q o criou
)
head(tx2gene_salmon)
#Isso tudo foi pra arrumar a diferença entre os nomes
txi <-tximport(arquivos_amostras,
               type = "salmon",
               tx2gene = tx2gene_salmon)
dim(txi$counts)
head(txi$counts)
#mostou que as colunas tao indo com o _Salmonado, vamo precisar tirar
colnames(txi$counts)      <- sub("_Salmonados$", "", colnames(txi$counts))
colnames(txi$abundance)   <- sub("_Salmonados$", "", colnames(txi$abundance))
colnames(txi$length)      <- sub("_Salmonados$", "", colnames(txi$length))
colnames(txi$counts)
#bonitinho sem o sufixo
#-------3 -> Data frames e coisarada
DadosAmostras <- data.frame(
  amostras = c("SRR21831518", "SRR21831519", "SRR21831520", 
               "SRR21837144", "SRR21837145", "SRR21837146"),
  tecido = c("cerebral", "cerebral", "cerebral",
             "adiposo", "adiposo", "adiposo")
)
row.names(DadosAmostras) <- DadosAmostras$amostras
#Data frame com as if que temos
#-----------4 -> Rodando o deseq2 
dds <- DESeqDataSetFromTximport(txi, colData = DadosAmostras, design =  ~ tecido)
dds
#--------5 -> Filtrando
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep]

#---botando a ref ent estou botando aqlq q eu espero que tenha menos expressao )
dds$tecido <- relevel(dds$tecido, ref = "cerebral")

dds <-DESeq(dds)
res <- results(dds)
summary(res)
#Fazendo com q o p value seja de 0.01 e n 0.1
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)
res0.05 <- results(dds,alpha = 0.05)
summary(res0.05)
