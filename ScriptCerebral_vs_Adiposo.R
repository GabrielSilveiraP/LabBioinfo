#Primeiro script tentando buscar a expressao diferencial
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(here)
library(ggrepel)
tabelagenes <- read.csv(here("prepDE", "tabelagenes.csv"), row.names = 1)
head(tabelagenes)
tabelatranscritos <- read.csv(here("prepDE", "tabelatranscritos.csv"), row.names =  1 )
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
ntd <-normTransform(dds)

res <- results(dds)

res

vsd <- vst(dds, blind = FALSE)
#Vst(Variance stabilazing Transformation)
#transformação considera o delineamento experimental (tecido cerebral vs adiposo)


#ja tamos nos resultados aqui

summary(res)
#Fazendo com q o p value seja de 0.01 e n 0.1
res0.01 <-  results(dds,alpha = 0.01)
res0.05 <-  results(dds,alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("id_completo") %>%
  mutate(
    id_anotacao = sub(".*\\|", "", id_completo),      # pega depois do | pipe
    id_stringtie = sub("\\|.*", "", id_completo),     # pega antes do |
    external_gene_name = id_anotacao                  # o nome já É a anotação, n precisa mais escrever algo
  )
summary(res0.01)
summary(res0.05)

sum(res$padj < 0.1, na.rm=TRUE)
#Plottando

plotMA(res0.01)
plotPCA(vsd, intgroup = "tecido") + 
  geom_text(aes(label = colnames(vsd)), vjust = 2, size = 3)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
#DadosAmostras <- as.data.frame(colData(ntd)[,c("tecido")])
pheatmap(assay(ntd)[select,],
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col = DadosAmostras[, "tecido", drop=FALSE])

#plotDispEsts(dds) #Esse sozinho é o fitted


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, "tecido", drop=FALSE])
nomes_linhas <- rownames(assay(ntd)[select,])
nomes_anotados <- res_df$external_gene_name[
  match(sub(".*\\|", "", nomes_linhas), res_df$id_anotacao)
]
labels <- ifelse(is.na(nomes_anotados), nomes_linhas, nomes_anotados)
pheatmap(assay(ntd)[select,],
         cluster_rows = FALSE,
         show_rownames = TRUE,   # <- aqui aparece o nome do gene
         labels_row = labels,    # <- aqui entra a anotação
         annotation_col = df)


# Ver quais são os 4 genes significativos
genes_significativos <- res_df %>% #Esse %>% é o pipe do tidyverse
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  select(id_completo, external_gene_name, log2FoldChange, padj)

head(genes_significativos)

select_sig <- which(rownames(assay(ntd)) %in% genes_significativos$id_completo)

# Labels com nome anotado
labels_sig <- sub(".*\\|", "", rownames(assay(ntd)[select_sig,]))

# Heatmap plottadinho
pheatmap(assay(ntd)[select_sig,],
         cluster_rows = FALSE, #É essa linha que altera a hierarquia, pelo menos a do lado
         show_rownames = TRUE,
         labels_row = labels_sig,
         annotation_col = DadosAmostras[, "tecido", drop=FALSE])

#fazendo outro df mas com os cabra bom de padj
res_df_ProVulcanPlot <- res_df[!is.na(res_df$padj), ]

#Oque ta abaixo é pra colorir os up e os down-coluna de cor
res_df_ProVulcanPlot$ExpressaoGenica <- "Não significativo"
res_df_ProVulcanPlot$ExpressaoGenica[res_df_ProVulcanPlot$padj < 0.05 & res_df_ProVulcanPlot$log2FoldChange > 1] <- "Up regulado (superexpresso)"
res_df_ProVulcanPlot$ExpressaoGenica[res_df_ProVulcanPlot$padj < 0.05 & res_df_ProVulcanPlot$log2FoldChange < -1] <- "Down regulado (subexpresso)"

#Fazendo label -> fazendo com que apareceça o nome dos melhores 5
top5genes <- head(res_df_ProVulcanPlot[order(res_df_ProVulcanPlot$padj, na.last = NA),"external_gene_name"], 5)
res_df_ProVulcanPlot$nomedosTop5 <- ifelse(res_df_ProVulcanPlot$external_gene_name %in% top5genes, res_df_ProVulcanPlot$external_gene_name, NA)

#Aq ta o core do volcano plot
ggplot(data = res_df_ProVulcanPlot, aes(x = log2FoldChange, y = -log10(pvalue), color = ExpressaoGenica, label = nomedosTop5)) + 
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed") +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values = c("red", "grey","blue"),
                     labels = c("Down regulado (subexpresso)","Não significativo","Up regulado (superexpresso)")) +
  geom_point()


