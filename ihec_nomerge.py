import pandas as pd
import sys
import re


def reorganize_cols(df_to_work, nci_file):

    df_final = get_ontology_nci(df_to_work, nci_file)
    # re-organizing columns
    cols = ['EpiRR', 'EpiRR_status', 'age', 'biomaterial_type', 'cell_type',
       'disease_ontology_term', 'donor_age_unit',
       'donor_id', 'donor_life_stage', 'health_state', 'line', 'markers',
       'project', 'sample_ontology_term', 'sex', 'taxon_id', 'tissue_type',
       'donor_health_status', 'donor_health_status_merge', 'donor_health_status_new',
       'disease', 'cause of death','donor_health_status_ontology_curie', 'Description_donor_health_status_ontology_curie',
       'donor_health_status_ontology_uri', 'Description_donor_health_status_ontology_uri' ,'disease_ontology',
       'disease_ontology_curie', 'Description_disease_ontology_curie',
       'disease ontology uri','disease_ontology_uri', 'Description_disease_ontology_uri']
    
    df_final = df_final[cols]
    # df_final.to_csv('IHEC_metadata_v07_modified_nomerge_version4.csv', index=False)


def get_ontology_nci(df_to_work, nci_file):

    dict_ncit, dict_ncim = create_dict_ncit_ncim(nci_file)

    list_health_curie = []
    list_health_uri = []
    list_disease_curie = []
    list_disease_uri = []

    for i,row in df_to_work.iterrows():
        
        #health_curie
        #health_curie = re.match(r'^\w+[\: | = ]+(C\d+)$',str(row['donor_health_status_ontology_curie']))
        health_curie = str(row['donor_health_status_ontology_curie']).split(':')[-1]
        map_dicts(health_curie, dict_ncit, dict_ncim, list_health_curie)
        
        #health_uri
        health_uri = str(row['donor_health_status_ontology_uri']).split('code=')[-1].split('&')[0] #some url with code=C123&type=restriction
        map_dicts(health_uri, dict_ncit, dict_ncim, list_health_uri)

        # disease_curie
        if 'http' in str(row['disease_ontology_curie']):
            disease_curie = str(row['disease_ontology_curie']).split('code=')[-1] #three cells
            map_dicts(disease_curie, dict_ncit, dict_ncim, list_disease_curie)

        else:
            disease_curie = str(row['disease_ontology_curie']).split(':')[-1]
            map_dicts(disease_curie, dict_ncit, dict_ncim, list_disease_curie)

        # disease_uri
        disease_uri = str(row['disease_ontology_uri']).split('code=')[-1].split('&')[0]
        map_dicts(disease_uri, dict_ncit, dict_ncim, list_disease_uri)
    
    df_to_work['Description_donor_health_status_ontology_curie'] = list_health_curie
    df_to_work['Description_donor_health_status_ontology_uri'] = list_health_uri
    df_to_work['Description_disease_ontology_curie'] = list_disease_curie
    df_to_work['Description_disease_ontology_uri'] = list_disease_uri

    return df_to_work

    
def map_dicts(var, dict_ncit, dict_ncim, list_to_app):
        
        if var in dict_ncit.keys():
            list_to_app.append(dict_ncit.get(var))
        elif var in dict_ncim.keys():
            list_to_app.append(dict_ncim.get(var))
        else:
            list_to_app.append('N/A')


def create_dict_ncit_ncim(nci_file):

    dict_ncit = dict(zip(nci_file['ncit'], nci_file['ncit_term'])) #123701
    dict_ncim = dict(zip(nci_file['ncim'], nci_file['ncim_term'])) #123695

    return dict_ncit, dict_ncim 



def map_ontologies(df_to_work):

    dict_ont = {'ncim:c0549184':'None', 'ncim:c115222':'Unknown', 'ncim:c0277545':'None',
                'ncit:c0549184':'None', 'ncit:c115222':'Unknown', 'ncit:c0277545':'None'}

    # print(df_to_work[df_to_work['disease'].str.contains('Blunt', na=False)])

    df_to_work['donor_health_status_new'] = df_to_work['donor_health_status_new'].fillna(df_to_work['donor_health_status_ontology_curie'].str.lower().map(dict_ont))
    df_to_work['donor_health_status_new'] = df_to_work['donor_health_status_new'].fillna(df_to_work['disease_ontology_curie'].str.lower().map(dict_ont))

    #IHECRE00003764.1
    # print('after',df_to_work[df_to_work['disease'].str.contains('Blunt', na=False)])
    return df_to_work


def to_map_health(df1,df2): #ok

    dict_to_map = dict_disease_cause_death(df1)
    # print(dict_to_map['IHECRE00003764.1'])
    # sys.exit()
    df_ori_col = get_desired_cols(df1) #v01 - correct
    df2.rename({'donor_health_status': 'donor_health_status_merge'}, axis=1, inplace=True) #v07
    df_to_work = df2.merge(df_ori_col, on='EpiRR' ,how='left') #merging ori columns
    df_to_work['donor_health_status_new'] = df_to_work['donor_health_status']
    df_to_work['donor_health_status_new'] = df_to_work['donor_health_status_new'].fillna(df_to_work['EpiRR'].map(dict_to_map))
    
    return df_to_work
    

def dict_disease_cause_death(df1): #to_map_health

    dict_to_map = {}

    for i,row in df1.iterrows():
        if pd.notna(row['disease']):
            dict_to_map[row['EpiRR']] = row['disease']

        elif pd.isna(row['disease']) and pd.notna(row['cause of death']):
            dict_to_map[row['EpiRR']] = row['cause of death']
    
    return dict_to_map


def get_desired_cols(df1):
    
    return df1[['EpiRR','donor_health_status','disease','cause of death','donor_health_status_ontology_curie',
    'donor_health_status_ontology_uri','disease_ontology','disease_ontology_curie','disease ontology uri',
    'disease_ontology_uri']]


def main():

    print('Starting...')
    df_raw = pd.read_csv(sys.argv[1]) #version 1 (EpiRR updated - no merged columns)
    df_v7 = pd.read_csv(sys.argv[2]) #version 7 (last version)
    nci_file = pd.read_csv(sys.argv[3], sep='|', index_col=False, names = ['ncit','ncim', 'ncit_term', 'ncim_term'])
    df_to_work = to_map_health(df_raw,df_v7) #correct
    df_to_work_ont = map_ontologies(df_to_work) #terms from curie columns
    df_to_work_ont.to_csv('frosi_ont.csv')
    sys.exit()
    
    reorganize_cols(df_to_work_ont, nci_file) #to create columns with term description and reorganize cols
    # get_ontology_nci(df_to_work_ont,nci_file) 
    
    
    
    print('Df saved!')
    
     

if __name__ == "__main__":



    main()