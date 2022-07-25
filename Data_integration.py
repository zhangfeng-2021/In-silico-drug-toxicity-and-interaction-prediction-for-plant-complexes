# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 22:00:35 2022

@author: User
"""
import pandas as pd
import time

start = time.time()
dftoxitotal = pd.DataFrame()
dfsimitotal = pd.DataFrame()
start = time.time() 
similarity = 1
interv = -0.11
while similarity <=1 and similarity >=0.10:
    print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
    with open(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
        dfsimilar = pd.read_csv(fin, index_col= 0) # Read similar compound dataframe
    with open(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Toxi_all_infor_similar_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
        dftoxi = pd.read_csv(fin, index_col= 0) # Read similar compound dataframe
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv", index_col= 0)
    dfinhouse = dfinhouse.append(dfsimilar, ignore_index = True)
    dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID"])
    dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv") 
    print("More similar compounds properties added to in-house vector")
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv", index_col= 0)
    dfinhouse = dfinhouse.append(dftoxi, ignore_index = True)
    dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID","Similarity_score"], keep='last')    
    dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv")
    print("More toxicity and drug interaction information added to in-house vector")
    dftoxitotal = dftoxitotal.append(dftoxi, ignore_index = True)
    dfsimitotal = dfsimitotal.append(dfsimilar, ignore_index = True)
    similarity += interv 
dftoxitotal = dftoxitotal.sort_values(['Active_compound_name', 'Similarity_score'], ascending=[True, False], axis = 0)
dfsimitotal = dfsimitotal.sort_values(['Active_compound_name', 'Similarity_score'], ascending=[True, False], axis = 0)
dftoxitotal = dftoxitotal.drop_duplicates(subset = ["Active_compound_name","PubChem_CID"], keep = "first")

dftoxitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv")
mask = (dftoxitotal["Drug Induced Liver Injury"].notna()|dftoxitotal["Hepatotoxicity"].notna() |dftoxitotal["Evidence for Carcinogenicity"].notna() |dftoxitotal["Carcinogen Classification"].notna())
dftoxitotal = dftoxitotal[mask]
print(dftoxitotal.shape)
dftoxitotal = dftoxitotal.sort_values(['Active_compound_name', 'Similarity_score'], ascending=[True, False], axis = 0)
dftoxitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum_similarity_threshold_analysis_data.csv")
dfsimitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\similar_comp_properties_sum.csv")
end = time.time()
print('execution time was {} minutes'.format(round((end-start)/60, 2)))       
        
        
        
