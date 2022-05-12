# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 21:50:53 2021

@author: User
"""
import pandas as pd

with open(r"Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)

similarity = 0.7 
RepreComList1= list(set(df1["Active compound name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active compound name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active compound name"]))
print(activenumber1)

similarity = 0.561 
RepreComList1= list(set(df1["Active compound name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active compound name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active compound name"]))
df3 = df2.drop_duplicates(subset=["Active compound name", "PubChem CID of similar compounds"], keep='first', inplace=False)
df3.to_csv(r"prediction_0.561.csv")
print(activenumber1)