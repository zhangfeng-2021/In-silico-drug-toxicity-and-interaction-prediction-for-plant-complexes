# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 09:46:09 2021

@author: User
"""
import pubchempy as pcp
import pandas as pd
import time,math
from requests_html import HTMLSession
from lxml.html import fromstring
import requests, re, os, argparse,json
import urllib.parse
import urllib.request
import numpy as np
import pyautogui as pag

similarity = 1
interv = -0.11
start1 = time.time()

while similarity <=1 and similarity >=0.10:
    df1 = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\similar_comp_pool_swiss_2022_07_08.csv")
    RepreComList = list(set(df1["Active_compound_name"]))[0:]
    # print(RepreComList)
    mask = (df1["Similarity_score"] > (similarity + interv)) & (df1["Similarity_score"] <= similarity) & (df1["Active_compound_name"].isin(RepreComList))
    df1 = df1[mask]
    # print(df1)
    compoundsmiles = df1['SMILES_of_similar']
    activenamelist = df1["Active_compound_name"]
    CHEMBLlist = df1["CHEMBL_Similar_Compound_ID"]
    Similarity_scorelist = df1["Similarity_score"]
    df2 = pd.DataFrame(columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])
    dfstore = pd.DataFrame(columns =['Similar not found in Pubchem',"SMILES_of_similar","CHEMBLID"])
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv")
    i = 0
    try:
        for smiles, name, CHEMBLID, prob in zip(compoundsmiles[0:], activenamelist[0:], CHEMBLlist[0:], Similarity_scorelist[0:]):
            a,b = pag.position()
            time2 = time.time()
            timecost1 = time2 - start1
            print(math.ceil(timecost1/30))
            if math.ceil(timecost1/30) %100 == 0:
                            time.sleep(20)
            else:
                print("=====================================")
                print(df1.shape)
                print(i+1, CHEMBLID)
                print(str((round((similarity+interv),3), round(similarity,3))))
            
            dfins = dfinhouse[dfinhouse["CHEMBL_Similar_Compound_ID"]== CHEMBLID]
            if CHEMBLID in dfinhouse["CHEMBL_Similar_Compound_ID"].values.tolist():
                print("This compound has been stored in former vector")
                # print(dfins)
                dfins = dfins.drop(labels = ["Active_compound_name","Similarity_score"], axis= 1)
                print("Drop two coloumns which are not useful")
                # print(dfins)
                dfins = dfins.drop_duplicates(subset= "CHEMBL_Similar_Compound_ID", keep='first')
                print("Drop duplicates")
                # print(dfins)
                dfins["Active_compound_name"] = name
                dfins["Similarity_score"] = prob
                df2 = df2.append(dfins, ignore_index = True)
                print("It has been stored successfully")
                i += 1
            else:
                try:
                    cidlist = pcp.get_cids(CHEMBLID,'name')
                    if cidlist and cidlist[0] != 0 :
                        # print("could by name")
                        # print(cidlist)
                        for cids in cidlist:
                            print("PubChem_CID: " + str(cids))
                            df = pd.DataFrame([[name, CHEMBLID, cids, prob]],columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])
                            df2 = df2.append(df, ignore_index = True)
                            print(time.time(), "Retrieve success by name")
                            i += 1
                    else:
                        cidlist = pcp.get_cids(smiles,'smiles')
                        if cidlist:
                            for cids in cidlist:
                                df = pd.DataFrame([[name, CHEMBLID, cids, prob]],columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])            
                                # print(str(cids))
                                df2 = df2.append(df, ignore_index = True)
                                print(time.time(),"Retrieve success by ID-",CHEMBLID)
                                i += 1
                        else:
                            print(name,smiles,'no hit', time.time())
                            dfstore.loc[i] = ([name, smiles, CHEMBLID])
                            i += 1
                except:
                    print(name,smiles,'no hit', time.time())
                    dfstore.loc[i] = ([name, smiles, CHEMBLID])
                    i += 1
                     
        retrinamelist = df2['Active_compound_name']
        df2.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv") 
        seta = set(activenamelist)-set(retrinamelist) 
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_compound_pubchem_not_found_"+str(round(similarity, 3))+".csv")    
        end1 = time.time()
        # print('Compound not found in Pubchem:', seta)
        print('execution time was {} minutes'.format(round((end1-start1)/60, 2)))
    except:
        retrinamelist = df2['Active_compound_name']
        df2.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv") 
        seta = set(activenamelist)-set(retrinamelist) 
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_compound_pubchem_not_found_"+str(round(similarity, 3))+".csv")    
        end1 = time.time()
        # print('Compound not found in Pubchem:', seta)
        print('execution time was {} minutes'.format(round((end1-start1)/60, 2)))
        raise
    similarity += interv  