# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:05:59 2022

@author: User
"""
import pandas as pd
import time

herbinteresting = "Spatholobus Suberectus Dunn"
with open(r"D:\Study at HKU small documents\Publication\in-silico toxicity and drug interaction prediction_ssp\submission\Briefings in bioinformatics\Raw data updated on 2022 07 15\Active compounds_TCMSP\Active_compounds_pool.csv") as fin:
    dfactive = pd.read_csv(fin, index_col = 0)
mask = (dfactive["Herb"]== herbinteresting)
dfactive = dfactive[mask]
activelist = dfactive["Molecule Name"]
print(dfactive.head())

with open(r"D:\Study at HKU small documents\Publication\in-silico toxicity and drug interaction prediction_ssp\submission\Briefings in bioinformatics\Results\Prediction result\Interaction Network\drug_interaction_0.6171.csv") as fin:
    dfinteraction = pd.read_csv(fin, index_col = 0)
mask = (dfinteraction["Active_compound_name"].isin(activelist))
dfinteraction = dfinteraction[mask]
dfinteraction.to_csv(r"D:\Study at HKU small documents\Publication\in-silico toxicity and drug interaction prediction_ssp\submission\Briefings in bioinformatics\Results\Prediction result\Interaction Network\drug_interaction_"+herbinteresting+".csv")
print(dfinteraction.head())
    