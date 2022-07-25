# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 21:50:53 2021

@author: User
"""
import pandas as pd

with open(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
    
similarity = 1
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))
    

similarity = 0.759 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.6171 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
df3 = df2.drop_duplicates(subset=["Active_compound_name", "PubChem_CID"], keep='first', inplace=False)
df3.to_csv(r"prediction_"+str(similarity)+".csv")
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))


with open(r"in-silico toxicity and drug interaction prediction\integration\similar_comp_properties_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
    
similarity = 1 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.759 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.6171 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

