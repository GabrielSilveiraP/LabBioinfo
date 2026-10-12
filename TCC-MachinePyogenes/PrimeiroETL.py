#%%
import pandas as pd

# %%
tabelona = pd.read_csv("Cópia de Formulário_Notificação_Streptococcus_pyogenes-VPURA.csv", sep=";", header= 1, encoding= "utf-8-sig")
# %%
print([c for c in ["Número interno",
                          "Número da amostra" ,
                          "Data resultado Laboratório Externo" ,
                          "Data envio Laboratório Externo" ,
                          "O evento trata-se de:",
                          "Leia as informações na descrição da seção acima:"	,
                          "Número da notificação SINAN:	Pontuação"] if c not in tabelona.columns])
#%%
#Tive que fazer isso porq essas n apagavam de jeito nenhum
palavras = ["Leia as informações", "Número da notificação SINAN","Carimbo de data/hora"]
cols_apagar = [c for c in tabelona.columns if any(p in c for p in palavras)]
print(cols_apagar)  # confira antes de apagar
tabelona = tabelona.drop(columns=cols_apagar)
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
                          "Nome completo do responsável pela notificação:"
                          ])
# %%
tabelona
# %%
preenchimento = pd.DataFrame({"preenchidos": tabelona.notna().sum(), "%": (tabelona.notna().mean() * 100).round(1)}).sort_values("%")
# %%
preenchimento
# %%
tabelona = tabelona.replace(
                {"Não": 0,
                "Sim": 1
})

#%%
tabelona
#%%
#Preciso arrumar a coluna de idade, que veio toda misturada, entre meses e anos. Noa criei 2 colunas, meses e anos por conta que temos pouquissimos casos de meses, dai o modelo via ler e ter q passar por muitas zeros desnecessariamente

limpo = tabelona["Idade (em anos):"].astype(str).str.strip().str.lower()

num   = pd.to_numeric(limpo.str.replace(",", "."), errors="coerce")        # só número: "11", "48"
anos  = limpo.str.extract(r"(\d+)\s*(?:anos?|a)\b")[0].astype(float)       # "57 anos", "8A", "6ANOS"
meses = limpo.str.extract(r"(\d+)\s*m")[0].astype(float)                   # "10 meses", "8m"
dias  = limpo.str.extract(r"(\d+)\s*d")[0].astype(float)                   # "10 dias"

tabelona["Idade (em anos):"] = num.fillna(anos.fillna(0) + meses.fillna(0)/12 + dias.fillna(0)/365)
#%%
#Arredondar p 2 casas decimais depois d virgula e tbm mostrar a contagem unica
tabelona["Idade (em anos):"] = tabelona["Idade (em anos):"].round(2)
tabelona["Idade (em anos):"].value_counts().sort_index()
#A coluna atualizacao evolucao (comumente evo2) todos os na irei preencher com a coluna anterior
tabelona["Atualização Evolução"] = tabelona["Atualização Evolução"].fillna(tabelona["Evolução do paciente:"])
#%%
tabelona["Atualização Evolução"] = tabelona["Atualização Evolução"].replace({
                "Óbito em 15/11": "Óbito",
                "Evadiu do hospital" : "Não se aplica",
                "Mantém internação em 28/11" : "Mantém internação",
                "Setor de internação" : "Alta hospitalar"
})
#%%
#Encontrei uma linha com uma informaçao que nao ta coerente, por isso comando abaixo do value.counts è uma pesquisa
tabelona["Atualização Evolução"].value_counts().sort_index()
#tabelona.loc[tabelona["Atualização Evolução"] == "Setor de internação"]
#Nao ta aqui porq botei o replace ali em cima
#%%
    #dummy_vars = [
    # "Se paciente internado, qual setor de internação?",
    # " Evolução do paciente:","Atualização Evolução"]

#%%
#df["tabelona"] = (df["Evolução do paciente:.óbito">0).astype(int)