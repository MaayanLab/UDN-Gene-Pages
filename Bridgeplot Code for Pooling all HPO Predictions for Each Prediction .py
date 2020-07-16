#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#! /usr/bin/python

import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from maayanlab_bioinformatics.plotting import bridge_plot

#read entire database into dataframe


#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']



#make dataframe for enrichr graph 
Ranks = []
Property = []
matching = []
c = 0

for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/enrichr/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    
    #store phenotypes for enrichr from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            #append for graph
            c+=1;
            Ranks.append(c)
            Property.append(dict['property'])
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    
    ##### GENERIF ####3##
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w]) 
        matching.append(val)
    
    #sort the list by ranks
    b.sort(key=lambda x:x[1])
    
    
    
    
    #ENRICHR GRAPH
    #make dataframe for graphs 
    graph = pd.DataFrame({'Property': a})
    
    #get only first 
    dff = graph.head(100)
    
    
    
    #make dataframe for matches
    a = pd.DataFrame({'HGO Disease ID': matching})


    # extract HP ids
    df_sorted = dff
    df_sorted["HP_ID"] = [x[1].replace(")", "") for x in df_sorted["Property"].str.split("(")]
    select = df_sorted['HP_ID']

    a = pd.DataFrame({'HGO Disease ID': matching})
    
    
    #
    #sort my data
    # my_data_sorted = my_data.sort_values('HPO Matches Using Autorif Co-Occurrence')
    # (my_data_sorted['HPO Matches Using Autorif Co-Occurrence'] == True).hist().plot()
    x, y = bridge_plot(df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    plt.plot(x, y,label="Gene:" + gene)

   


    
    #add labels
    plt.xlabel('Ranked Phenotype Predictions from Geneshot')
    plt.ylabel('# of Hits')
    plt.title('UDN GENE: ' + gene + '  Matches for HPO Terms From GeneShot Prediction Methods')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.vlines(np.where((df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    , 1, 0), ymin=-1, ymax=0)
    
    

plt.show()
    
   


# In[13]:


#! /usr/bin/python

import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from maayanlab_bioinformatics.plotting import bridge_plot

#read entire database into dataframe


#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']



#make dataframe for generif graph 
generif_Ranks = []
generif_Property = []
generif_matching = []
generif_c = 0





for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    
    
    # GENERIF DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    GENERIF_BACK_URL = '/library/Human_Phenotype_Ontology/similarity/generif/offset/0/limit/200'
    GENERIF_URL = URL + gene + GENERIF_BACK_URL 

    GENERIF_response = requests.get(
    GENERIF_URL
 )
    if not GENERIF_response.ok:
        raise Exception('Error during query')

    GENERIF_data = json.loads(GENERIF_response.text)


    GENERIF_results = GENERIF_data.get('results')
    
    if GENERIF_results is None:
        continue

    

    
    #store phenotypes from geneshot data
    generif_a = []
    generif_rank = 1
    
    for dict in GENERIF_results:
        generif_rank+=1
        for key,value in dict.items():
            #append for graph
            generif_c+=1;
            generif_Ranks.append(generif_c)
            generif_Property.append(dict['property'])
            generif_a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    generif_x = set(phenotype) & set(generif_a)
    

    
    
    #create a list of matches with ranks for generif data
    generif_b = []
    generif_w = 0;
    for val in generif_x:
        generif_w = generif_a.index(val);
        generif_w+=1;
        z= 'Rank:' + str(w)
        
        #these will used for creating new dataframe
        generif_b.append([val, generif_w]) 
        generif_matching.append(val)
    
    #sort the list by ranks
    generif_b.sort(key=lambda x:x[1])
    

 
    
    #GENERIF GRAPH
    #make dataframe for graphs
    generif_graph = pd.DataFrame({'Property': generif_a})
    
    #get only first 
    generif_dff = generif_graph.head(100)
    
    
    
    #make dataframe for matches
    generif_ab= pd.DataFrame({'HGO Disease ID': generif_matching})
    
   


    # extract HP ids
    generif_df_sorted = generif_dff 
    generif_df_sorted["HP_ID"] = [x[1].replace(")", "") for x in generif_df_sorted["Property"].str.split("(")]
    generif_select = generif_df_sorted['HP_ID']
    
   
    
    #generif graph
    x2, y2 = bridge_plot(generif_df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in generif_ab["HGO Disease ID"].dropna().astype(str).str.split("(")]))
    plt.plot(x2, y2, label="Generif")
   


    
    #add labels
    plt.xlabel('Ranked Phenotype Predictions from Geneshot')
    plt.ylabel('# of Hits')
    plt.title('UDN GENE Matches for HPO Terms From GeneShot Generif Methods')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.vlines(np.where((df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    , 1, 0), ymin=-1, ymax=0)
    
    
   


# In[21]:



import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from maayanlab_bioinformatics.plotting import bridge_plot

#read entire database into dataframe


#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']



#make dataframe for generif graph 
autorif_Ranks = []
autorif_Property = []
autorif_matching = []
autorif_c = 0

#make dataframe for generif graph 
tagger_Ranks = []
tagger_Property = []
tagger_matching = []
tagger_c = 0



for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/enrichr/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    
    
    # AUTORIF_ DATA
    AUTORIF_URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    AUTORIF_BACK_URL = '/library/Human_Phenotype_Ontology/similarity/autorif/offset/0/limit/200'
    AUTORIF_URL = AUTORIF_URL + gene + AUTORIF_BACK_URL 

    AUTORIF_response = requests.get(
        AUTORIF_URL
 )
    if not AUTORIF_response.ok:
        raise Exception('Error during query')

    AUTORIF_data = json.loads(AUTORIF_response.text)


    AUTORIF_results = AUTORIF_data.get('results')
    
    if AUTORIF_results is None:
        continue
        
    
     
        
     #AUTORIF DATA 
    
    #store phenotypes from geneshot data
    autorif_a = []
    autorif_rank = 1
    
    for dict in AUTORIF_results:
        autorif_rank+=1
        for key,value in dict.items():
            #append for graph
            autorif_c+=1;
            autorif_Ranks.append(autorif_c)
            autorif_Property.append(dict['property'])
            autorif_a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    #phenotype is my UDN table phenotypes
    autorif_x = set(phenotype) & set(autorif_a)
    

    
    
    #create a list of matches with ranks for autorif n data
    autorif_b = []
    autorif_w = 0;
    for val in autorif_x:
        autorif_w = autorif_a.index(val);

        autorif_w+=1;
        z= 'Rank:' + str(w)
        
        #these will used for creating new dataframe
        autorif_b.append([val, autorif_w]) 
        autorif_matching.append(val)
    
    #sort the list by ranks
    autorif_b.sort(key=lambda x:x[1])
    
    
    
    
    

    
    #AUTORIF GRAPH
    #make dataframe for graphs
    autorif_graph = pd.DataFrame({'Property': autorif_a})
    
    #get only first 
    autorif_dff = autorif_graph.head(75)
    
    
    
    #make dataframe for matches
    autorif_ab= pd.DataFrame({'HGO Disease ID': autorif_matching})
    
   
    # extract HP ids
    autorif_df_sorted = autorif_dff
    autorif_df_sorted["HP_ID"] = [x[1].replace(")", "") for x in autorif_df_sorted["Property"].str.split("(")]
    autorif_select = autorif_df_sorted['HP_ID']
    
    

    
    #autorif graph
    x3, y3 = bridge_plot(autorif_df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in autorif_ab["HGO Disease ID"].dropna().astype(str).str.split("(")]))
    plt.plot(x3, y3, label="Autorif")
   

    
    #add labels
    plt.xlabel('Ranked Phenotype Predictions from Geneshot')
    plt.ylabel('# of Hits')
    plt.title('Geneshot Autorif HPO Matches for UDN Genes')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.vlines(np.where((df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    , 1, 0), ymin=-1, ymax=0)


# In[24]:



import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from maayanlab_bioinformatics.plotting import bridge_plot

#read entire database into dataframe


#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']




#make dataframe for generif graph 
tagger_Ranks = []
tagger_Property = []
tagger_matching = []
tagger_c = 0



for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    
      
    # TAGGER DATA
    TAGGER_URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    TAGGER_BACK_URL = '/library/Human_Phenotype_Ontology/similarity/tagger/offset/0/limit/200'
    TAGGER_URL = TAGGER_URL + gene + TAGGER_BACK_URL 

    TAGGER_response = requests.get(
        TAGGER_URL
 )
    if not TAGGER_response.ok:
        raise Exception('Error during query')

    TAGGER_data = json.loads(TAGGER_response.text)


    TAGGER_results = TAGGER_data.get('results')
    
    if TAGGER_results is None:
        continue

    
   
    
    #TAGGER
    #store phenotypes from geneshot data
    tagger_a = []
    tagger_rank = 1
    
    for dict in TAGGER_results:
        tagger_rank+=1
        for key,value in dict.items():
            #append for graph
            tagger_c+=1;
            tagger_Ranks.append(autorif_c)
            tagger_Property.append(dict['property'])
            tagger_a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    #phenotype is my UDN table phenotypes
    tagger_x = set(phenotype) & set(tagger_a)
    

    
    
    #create a list of matches with ranks for autorif n data
    tagger_b = []
    tagger_w = 0;
    for val in tagger_x:
        tagger_w = tagger_a.index(val);
        tagger_w+=1;
        z= 'Rank:' + str(w)
        
        #these will used for creating new dataframe
        tagger_b.append([val, tagger_w]) 
        tagger_matching.append(val)
    
    #sort the list by ranks
    tagger_b.sort(key=lambda x:x[1])
    
    
    
    

    #TAGGER
    
    #make dataframe for graphs
    tagger_graph = pd.DataFrame({'Ranks': tagger_Ranks, 'Property': tagger_Property})
    
    #get only first 
    tagger_dff = tagger_graph.head(75)
    
    
    
    #make dataframe for matches
    tagger_ab= pd.DataFrame({'HGO Disease ID': tagger_matching})
    
   
    # extract HP ids
    tagger_df_sorted = tagger_dff
    tagger_df_sorted["HP_ID"] = [x[1].replace(")", "") for x in tagger_df_sorted["Property"].str.split("(")]
    tagger_select = tagger_df_sorted['HP_ID']
    
    
   
    #tagger graph
    x4, y4 = bridge_plot(tagger_df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in tagger_ab["HGO Disease ID"].dropna().astype(str).str.split("(")]))
    plt.plot(x4, y4, label="Tagger")
   


    #add labels
    plt.xlabel('Ranked Phenotype Predictions from Geneshot')
    plt.ylabel('# of Hits')
    plt.title('UDN GENE: ' + gene + '  Matches for HPO Terms From GeneShot Prediction Methods')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.vlines(np.where((df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    , 1, 0), ymin=-1, ymax=0)
    


# In[28]:


import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from maayanlab_bioinformatics.plotting import bridge_plot

#read entire database into dataframe


#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']




#make dataframe for archs4 graph 
archs4_Ranks = []
archs4_Property = []
archs4_matching = []
archs4_c = 0




for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
   
     #ARCHS4
    # ARCHS4 DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/coexpression/offset/0/limit/200'
    ARCHS4_URL = URL + gene + BACK_URL
    
    ARCHS4_response = requests.get(
    ARCHS4_URL
 )
    if not ARCHS4_response.ok:
        raise Exception('Error during query')
    
    ARCHS4_data = json.loads(ARCHS4_response.text)
    
    ARCHS4_results = ARCHS4_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ARCHS4_results is None:
        continue
    

    
    
    #ARCHS4
    #store phenotypes from geneshot data
    archs4_a = []
    archs4_rank = 1
    
    for dict in ARCHS4_results:
        tagger_rank+=1
        for key,value in dict.items():
            #append for graph
            archs4_c+=1;
            archs4_Ranks.append(archs4_c)
            archs4_Property.append(dict['property'])
            archs4_a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    #phenotype is my UDN table phenotypes
    archs4_x = set(phenotype) & set(archs4_a)
    

    
    
    #create a list of matches with ranks for autorif n data
    archs4_b = []
    archs4_w = 0;
    for val in archs4_x:
        archs4_w = archs4_a.index(val);
        archs4_w+=1;
        z= 'Rank:' + str(w)
        
        #these will used for creating new dataframe
        archs4_b.append([val, archs4_w]) 
        archs4_matching.append(val)
    
    #sort the list by ranks
    archs4_b.sort(key=lambda x:x[1])
    
    
    
    #ARCHS4
    #make dataframe for graphs
    archs4_graph = pd.DataFrame({'Property': archs4_a})
    
    #get only first 
    archs4_dff = archs4_graph.head(100)
    
    
    
    #make dataframe for matches
    archs4_ab= pd.DataFrame({'HGO Disease ID': archs4_matching})
    
   
    # extract HP ids
    archs4_df_sorted = archs4_dff
    archs4_df_sorted["HP_ID"] = [x[1].replace(")", "") for x in archs4_df_sorted["Property"].str.split("(")]
    archs4_select = archs4_df_sorted['HP_ID']
    
    
    
    

   
    #archs4
    x5, y5 = bridge_plot(archs4_df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in archs4_ab["HGO Disease ID"].dropna().astype(str).str.split("(")]))
    plt.plot(x5, y5, label="ARCHS4")


    
    #add labels
    plt.xlabel('Ranked Phenotype Predictions from Geneshot')
    plt.ylabel('# of Hits')
    plt.title('UDN GENE: ' + gene + '  Matches for HPO Terms From GeneShot Prediction Methods')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.vlines(np.where((df_sorted["HP_ID"].isin([x[1].replace(")", "") for x in a["HGO Disease ID"].dropna().str.split("(")]))
    , 1, 0), ymin=-1, ymax=0)
    


# In[ ]:




