#Primeiro script tentando buscar a expressao diferencial
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)
listMarts()
#renv junto da here, setwd = um cmd p facilitar o caminho
tabelagenes <- read.csv("/home/gabriel/Documents/CoisasdoLab/QuantificacaoStringTIE-TabelasprepDE/tabelagenes.csv", row.names = 1 )
head(tabelagenes)
tabelatranscritos <- read.csv("/home/gabriel/Documents/CoisasdoLab/QuantificacaoStringTIE-TabelasprepDE/tabelatranscritos.csv", row.names =  1 )
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
  countData = tabelagenes,
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
vsd <- vst(dds, blind = FALSE)
ntd <- normTransform(dds)
DistanciaAmostras <- dist(t(assay(vsd)))

res


#ja tamos nos resultados aqui

summary(res)
#Fazendo com q o p value seja de 0.01 e n 0.1
res0.01 <-  results(dds,alpha = 0.01)
res0.05 <- results(dds,alpha = 0.05)


summary(res0.01)
summary(res0.05)

#Mais para frente foi notado que o padj estava muito alto, a linha abaixo corrige isso
resPadjAjustado <- subset(res0.05, padj < 0.05)
resPadjAjustado <- resPadjAjustado [order(resPadjAjustado$padj),]
head(resPadjAjustado)
#Plottando
plotDispEsts(dds)
plotMA(res) 
plotMA(res0.01)

SampleDistMatrix <- as.matrix(DistanciaAmostras)
rownames(SampleDistMatrix) <- paste(vsd$tecido, colnames(assay(vsd)), sep = " - ")
colnames(SampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(SampleDistMatrix,
         clustering_distance_rows=DistanciaAmostras,
         clustering_distance_cols=DistanciaAmostras,
         col=colors)

select <- order(rowMeans(counts(dds, normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, "tecido", drop=FALSE])
pheatmap(assay(ntd)[select,],
         cluster_rows = FALSE,
         show_rownames = TRUE,
         annotation_col = df)

# Ver se esse gene é estatisticamente significativo no DESeq2
         
res[grep("TcYC6_0172900", rownames(res)), ]

#depois de tudo isso, tivemos 4 genes que tiveram as melhores coberturas
genes_maisExpressos <-rownames(resPadjAjustado)
pheatmap(assay(ntd)[genes_maisExpressos,],
         cluster_rows = FALSE,
         show_rownames = TRUE,
         annotation_col = df)
#blast contra os transcritos anotados da ref
#blastn ( nucleiotideo contra nucleotideo)
#vulcan plot pelo valor de confiança
#Ortofinder
#usar os dados de cultura que o eric usa como base


library(ggplot2)

# Transformar em data.frame
res_df_volcano <- as.data.frame(res0.05) |>
  rownames_to_column("gene") |>
  filter(!is.na(padj))  # remover NAs

# Classificar os genes
res_df_volcano <- res_df_volcano |>
  mutate(significativo = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

# Plot
ggplot(res_df_volcano, aes(x = log2FoldChange, y = -log10(padj), color = significativo)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(data = filter(res_df_volcano, significativo != "NS"),
            aes(label = gene), vjust = -0.5, size = 3) +
  labs(
    title = "Volcano Plot - Cerebral vs Adiposo",
    x = "log2 Fold Change",
    y = "-log10(padj)",
    color = "Expressão"
  ) +
  theme_classic()

