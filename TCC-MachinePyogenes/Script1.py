# %%
import pandas as pd
from pysus import sih
#Importando oque  precisa

# %%
df_SP = sih(state="SP", year=2018, month=list(range(1,6)))

df_SC = sih(state="SC", year=2018, month=list(range(1,7)))
df_Concatenado = pd.concat([df_SP,df_SC])
#Baixando dados

# %%
df_Concatenado.head()
# %%
df_Concatenado.columns.tolist()
# %%
df_Metricas = df_Concatenado[[
    "SEXO",
    "RACA_COR",
    "ETNIA",
    "QT_DIARIAS",
    "DIAS_PERM",
    "MARCA_UTI",
    "CAR_INT",
    "COMPLEX",
    "UTI_MES_TO",
    "UTI_INT_TO",
    "DIAG_PRINC",
    "DIAG_SECUN",
    "MORTE",
    "CID_NOTIF",
    "DIAGSEC1",
    "DIAGSEC2",
    "CID_ASSO",
    "CID_MORTE",
    "PROC_REA",
    "IDENT",
    "SEQ_AIH5"
]]
#%%
df_Metricas.info()
#Aq eu vi que todos estavam como object, oque ta me impedindo de mexer como quiser
# %%
df_Metricas["ETNIA"] > 0
# Nao funciona porq ta em object
#%%
(df_Metricas["ETNIA"] == "0").sum()
#Só para checar quantos 0000 tem na lista
#%%
df_Metricas["ETNIA"].isna().sum()
df_Metricas["ETNIA"].dtype
#%%
df_Metricas["ETNIA"].unique()
#Pra ver os numeros diferentes - unicos- que tem nessa coluna, tem muito mais do que eu achava e tem erros de digitação do 0.
#%%
#Transformador0000 = 0000 isso aq substitui oq ta dentro dos (). Mas oque foi feito aqui é a transofrmação de tudo que era na, leia-se NaN, em 000 que á a medida que utilizamos quando não é falado a etnia do paciente.
df_Metricas.loc[:, "ETNIA"] = df_Metricas["ETNIA"].fillna("0000")

#%%
df_Metricas["DIAS_PERM"].isna().sum()
#Vi que tinha uma cacetada de Nan, entao rodei para que todas fossem substituidas por 0. Fiz isso com QT_DIARIAS, DIAS_PERM, UTI_MES_TO e UTI_INT_TO
#%%
df_Metricas.loc[:,["DIAS_PERM", "UTI_MES_TO","UTI_INT_TO" ]] = (df_Metricas[["DIAS_PERM", "UTI_MES_TO","UTI_INT_TO"]].fillna(0))
#%%
df_Metricas["UTI_MES_TO"] = df_Metricas["UTI_MES_TO"].fillna(0)
#%%
df_Metricas["UTI_INT_TO"] = df_Metricas["UTI_INT_TO"].fillna(0)
#%%
df_Metricas.loc[:, ["DIAS_PERM", "UTI_MES_TO", "UTI_INT_TO", "QT_DIARIAS"]] = (df_Metricas[["DIAS_PERM", "UTI_MES_TO", "UTI_INT_TO", "QT_DIARIAS"]].replace({None:0}).astype(int))
#Transformei esses em int, apra que eu consiga fazer contas. Nota: podia ter feito dessa maneira concatenada pra ali p cima.
#%%
df_Metricas["DIAS_PERM"].apply(type).value_counts()
#Retorna qual tipo, se ta com coisa misturada e value.count é pra contar essa qtd de vezes
#%%
for col in ["UTI_MES_TO", "UTI_INT_TO", "QT_DIARIAS"]:
    print(col, df_Metricas[col].apply(type).value_counts())
#Pra ver se realmente todos estavam com os numeros certos, cmo int e porque nao estava dando para fazer tudo junto
#%%
for col in ["DIAS_PERM", "UTI_MES_TO", "UTI_INT_TO", "QT_DIARIAS"]:
    df_Metricas[col] = df_Metricas[col].astype(int)
#%%
df_Metricas
#%%
#A primeira vez que fiz faltou colocar o & entre os filtros, lembrar de colocar ele ou | 
cids = ["J02.0", "J03.0", "A40.0", "B95.0"]
df_FiltradoCids = df_Concatenado[
    (df_Concatenado["DIAG_PRINC"].isin(cids)) |
    (df_Concatenado["DIAG_SECUN"].isin(cids)) |                   
    (df_Concatenado["DIAGSEC1"].isin(cids)) 
]

# %%
print(df_FiltradoCids.head())