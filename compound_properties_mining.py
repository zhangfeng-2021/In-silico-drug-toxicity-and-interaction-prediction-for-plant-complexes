# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 21:30:46 2021

@author: User
"""
import pubchempy as pcp
import pandas as pd
import time


start = time.time()
with open(r"in-silico toxicity and drug interaction prediction\Active compounds_TCMSP\Active_compounds_2022_07_08.csv") as fin:
    df1 = pd.read_csv(fin)
compoundname= df1['Molecule Name']
df2 = pd.DataFrame()

try: 
    
    for com in compoundname:
        c = pcp.get_compounds(com,'name')
        for i in c:
                df3 = pcp.compounds_to_frame(i)
                df3['Molecule Name'] = com
                df2 = df2.append(df3)
                print(com,i.isomeric_smiles)
                
except:
            print(com,'no hit')
            df2 = df2.append([com, 'no hit'])

retrinamelist = df2['Molecule Name']
df2.to_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties2022_07_08.csv") 
seta = set(compoundname)-set(retrinamelist) 
df4 = pd.DataFrame(list(seta), columns=['Compound not found in Pubchem'])    
df4.to_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties_no_retrie2022_07_08.csv")    
end = time.time()

print('Compound not found in Pubchem:', seta)
print('execution time was {} minutes'.format(round((end-start)/60, 2)))
