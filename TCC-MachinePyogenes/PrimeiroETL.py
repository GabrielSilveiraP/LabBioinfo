#%%
import pandas as pd
# %%
tabelona = pd.read_csv("Puras/Cópia de Formulário_Notificação_Streptococcus_pyogenes-VPURA.csv", sep=";", header= 1, encoding= "utf-8-sig")
# %%
#Limpando a planilha, remove linhas e colunas 100% vazias
tabelona = tabelona.dropna(how="all").dropna(axis=1, how="all")  
tabelona = tabelona.drop(columns=["Número interno",
                          "Número da amostra" ,
                          "Data resultado Laboratório Externo" ,
                          "Data envio Laboratório Externo" ,
                          "O evento trata-se de:",
                          "Pontuação",
                          "Endereço de e-mail",
                          "Nome completo do responsável pela notificação:",                          
                          "Nome da Instituição:",                          
                          "Amostras a serem enviadas para sequenciamento",                          
                          "Se outro, informe nome e CNES da unidade:"
                          ])
tabelona["Atualização Evolução"] = tabelona["Atualização Evolução"].replace({
                "Óbito em 15/11": "Óbito",
                "Evadiu do hospital" : "Não se aplica",
                "Mantém internação em 28/11" : "Mantém internação",
                "Setor de internação" : "Alta hospitalar"
})
#%%
palavras = ["Leia as informações", "Número da notificação SINAN", "Carimbo de data/hora"]
cols_apagar = [c for c in tabelona.columns if any(p in c for p in palavras)]
print(cols_apagar) 
tabelona = tabelona.drop(columns=cols_apagar)
#Fiz com que sempre que tivesse "Nao se aplica" na coluna de Atualizaçao evolucao fosse sobrescrito pela ultima categoria dada
#A coluna atualizacao evolucao (comumente evo2) todos os na irei preencher com a coluna anterior
tabelona["Atualização Evolução"] = (
    tabelona["Atualização Evolução"]
    .replace("Não se aplica", pd.NA)
    .fillna(tabelona["Evolução do paciente:"])
)
#%%
#Preciso arrumar a coluna de idade, que veio toda misturada, entre meses e anos. Nao criei 2 colunas, meses e anos por conta que temos pouquissimos casos de meses, dai o modelo via ler e ter q passar por muitas zeros desnecessariamente
col = tabelona["Idade (em anos):"]
limpo = col.astype(str).str.strip().str.lower()

num = pd.to_numeric(limpo.str.replace(",", ".", regex=False), errors="coerce")

anos = limpo.str.extract(r"(\d+)\s*(?:anos?|a)\b")[0].astype(float)
meses = limpo.str.extract(r"(\d+)\s*m")[0].astype(float)
dias = limpo.str.extract(r"(\d+)\s*d")[0].astype(float)

convertido = anos.fillna(0) + meses.fillna(0)/12 + dias.fillna(0)/365
tabelona["Idade (em anos):"] = num.fillna(convertido).round(2)
# %%
# tabelona['Carimbo de data/hora'] = pd.to_datetime(
#     tabelona['Carimbo de data/hora'], format='%m/%d/%Y %H:%M:%S'
# )
# print(tabelona['Carimbo de data/hora'].isna().sum())   # esperado: 0
# print(tabelona['Carimbo de data/hora'].min(), tabelona['Carimbo de data/hora'].max())

# %%
print(tabelona.duplicated().sum(), "linhas totalmente duplicadas")
print(tabelona["Atualização Evolução"].value_counts(dropna=False))
print(pd.crosstab(tabelona["Evolução do paciente:"],
                  tabelona["Atualização Evolução"], dropna=False))

# %%
tabelona.to_excel("TabelaLimpaClaudeEtapa1.xlsx", index=False)

# %%
