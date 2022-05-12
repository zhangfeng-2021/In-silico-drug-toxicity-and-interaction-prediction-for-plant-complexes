# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 11:58:08 2021

@author: User
"""
import pandas as pd
import ast
import os

df2 = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity score",
                              "Interactions","Drug_interacted", "Activities","Roles"])
df3 = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity score",
                              "Interactions","Drug_interacted", "Activities","Roles"])

with open(r"distinct drug_interaction_0.561.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
    container = [] 
for inter, cid, name, score in zip(df1["Interactions"], df1["PubChem CID of similar compounds"], 
                                   df1["Active_compound_name"], df1["Similarity score"]):
    inter = ast.literal_eval(inter)                                                            
    for j in inter:
        print(j)
        container.append(name)
        container.append(cid)
        container.append(name)        
        container.append(score)
        container.append(j)
        df3 = pd.DataFrame([container], columns = ["Active_compound_name","PubChem CID of similar compounds",
                                                   "Similarity score","Interactions"])
        df2 = df2.append(df3, ignore_index= True)
        container = [] 
df2.to_csv(r"drug_interaction_network_0.561.csv")