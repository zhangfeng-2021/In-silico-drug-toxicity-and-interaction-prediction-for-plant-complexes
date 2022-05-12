# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 21:30:46 2021

@author: User
"""
import pubchempy as pcp
import pandas as pd
import time


start = time.time()
df1 = pd.read_csv()
compoundsmiles = df1['SMILES_of_similar']
activenamelist = df1["Active_compound_name"]
CHEMBLlist = df1["CHEMBL_Similar_Compound_ID"]
scorelist = df1["Similarity_score"]
df2 = pd.DataFrame(columns=['Active_compound_name'])
df4 = pd.DataFrame(columns=['Similar not found in Pubchem',"SMILES_of_similar"])
i = 0
for smiles, name, CHEMBLID, score in zip(compoundsmiles[0:], activenamelist[0:], CHEMBLlist[0:], scorelist[0:]): 
        try:
            c = pcp.get_compounds(smiles,'smiles')
            for com in c:
                print(com)
                df3 = pcp.compounds_to_frame(com)
                df3['Active_compound_name'] = name
                df3["Similarity_score"] = score
                df2 = df2.append(df3)
                print(name,com.isomeric_smiles)
        except: 
            c = pcp.get_compounds(CHEMBLlist,'name')
            if len(c):
                for com in c:
                    print(com)
                    df3 = pcp.compounds_to_frame(com)
                    df3['Active_compound_name'] = name
                    df3["Similarity_score"] = score
                    df2 = df2.append(df3)
                    print(name,com.isomeric_smiles)
            else:
                print(name,smiles,'no hit')
                df4.loc[i] = ([name, smiles])
                i += 1
    
retrinamelist = df2['Active_compound_name']
df2.to_csv() 
seta = set(activenamelist)-set(retrinamelist) 

df4.to_csv()    
end = time.time()

print('Compound not found in Pubchem:', seta)
print('execution time was {} minutes'.format(round((end-start)/60, 2)))
