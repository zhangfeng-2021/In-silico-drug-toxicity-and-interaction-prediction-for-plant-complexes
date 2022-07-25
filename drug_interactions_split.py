# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 11:58:08 2021

@author: User
"""
import pandas as pd
import ast
import os

threshold = 0.561
dfintera = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity_score",
                              "Interactions","Drug_interacted", "Activities","Roles"])
dfins = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity_score",
                              "Interactions","Drug_interacted", "Activities","Roles"])

with open(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
container = []
df1 = df1[df1["Interactions"].notna()]
mask = df1["Similarity_score"] > threshold
df1 = df1[mask] 
for inter, cid, name, score in zip(df1["Interactions"], df1["PubChem_CID"], 
                                   df1["Active_compound_name"], df1["Similarity_score"]):
    inter = ast.literal_eval(inter)                                                            
    for j in inter:
        print(j)
        container.append(name)
        container.append(cid)
        container.append(score)
        container.append(j)
        dfins = pd.DataFrame([container], columns = ["Active_compound_name","PubChem CID of similar compounds",
                                                   "Similarity_score","Interactions"])
        dfintera = dfintera.append(dfins, ignore_index= True)
        container = []
dfintera.drop_duplicates(subset=(["Active_compound_name", "PubChem CID of similar compounds","Interactions"]), keep= 'first')
dfintera.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Drug_interaction\drug_interaction_"+str(threshold)+".csv")
print("Interaction data splitted")