# %%
import pandas as pd
from pysus import sih

#Importando oque  precisa

# %%
df_SP = sih(state="SP", year=2018, month=list(range(1,7)))

df_SC = sih(state="SC", year=2018, month=list(range(1,7)))
df_Concatenado = pd.concat([df_SP,df_SC])
#Baixando dados

# %%
print(df_Concatenado.head())
# %%
df_Concatenado[["UF_ZI", "MUNIC_RES", "IDADE", "SEXO", "RACA_COR", "ETNIA", "QT_DIARIAS", "MARCA_UTI", "UTI_MES_AL", "UTI_MES_TO", "UTI_INT_TO", "DIAG_PRINC", "DIAG_SECUN", "MORTE", "CID_NOTIF", "DIAGSEC1", "DIAGSEC2", "CID_ASSO", "CID_MORTE", "CONTRACEP1", "CONTRACEP2"]]
# %%
#A primeira vez que fiz faltou colocar o & entre os filtros, lembrar de colocar ele ou | 
cids = ["J02.0", "J03.0", "A40.0", "B95.0"]
df_FiltradoCids = df_Concatenado[
    (df_Concatenado["DIAG_PRINC"].isin(cids)) |
    (df_Concatenado["DIAG_SECUN"].isin(cids)) |                   
    (df_Concatenado["DIAGSEC1"].isin(cids)) 
]

# %%
print(df_FiltradoCids.head())