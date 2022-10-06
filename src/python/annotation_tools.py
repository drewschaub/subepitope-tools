def fetch_jsons(id_list, json_path, target='antigen'):
    """Fetch JSONs

    Args:
        id_list: most likely 4 character PDB ID codes
        target: string describing target (e.g. spike protein, fusion peptide, gp41)
        json_path: where json files will be stored

    Yields:
        df: dataframe object describing jsons
    """
    json_filenames = []
    error_list = []
    for i in tqdm(range(PDB_count)):
        PDB = id_list[i]
        request_url = api_url + PDB
        print(request_url)
        response = requests.get(request_url)

        polymer_entity_ids = response.json()["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
        
        for polymer_entity_id in polymer_entity_ids:
            polymer_entity_url = "https://data.rcsb.org/rest/v1/core/polymer_entity/" + PDB + "/" + polymer_entity_id
            
            # Get features of interest as json, and save json to avoid future scraping
            json_filename = "{}_{}_{}_polymer_entitiy.json".format(PDB, polymer_entity_id, target)
            json_full_path = Path(json_path, json_filename)
            json_filenames.append(json_filename)
            
            response = requests.get(polymer_entity_url)
            Path(json_path, json_filename).write_bytes(response.content)
            
            # Initialize a dictionary to store our features of interest
            polymer_entity_dict = {'PDB' : PDB,
                                'polymer_entity_id' : polymer_entity_id, 
                                'URL' : polymer_entity_url,
                                'json_filename' : json_filename
                                }

            polymer_entity_dict['pdbx_description'] = response.json()['rcsb_polymer_entity']['pdbx_description']
            polymer_entity_dict['provenance_code'] = response.json()['rcsb_polymer_entity']['rcsb_macromolecular_names_combined'][0]['provenance_code']
            
            json = response.json()['entity_poly']
            keys = ['pdbx_seq_one_letter_code', 'pdbx_strand_id', 'rcsb_sample_sequence_length']
            values = get_json_features(json, keys)
            polymer_entity_dict.update(dict(zip(keys, values)))
            polymer_entity_dict['pdbx_strand_id_count'] = polymer_entity_dict['pdbx_strand_id'].count(',') + 1
            
            keys = ['pdbx_gene_src_scientific_name']
            try:
                json = response.json()['entity_src_gen'][0]
                values = get_json_features(json, keys)
                # Some scientific names come back all capitalized. This will standardize them
                values = [value.capitalize() for value in values]
                polymer_entity_dict.update(dict(zip(keys, values)))
            except:
                values = [None for key in keys]
                polymer_entity_dict.update(dict(zip(keys, values)))
            
            keys = ['feature_id', 'name', 'provenance_source', 'type']
            try:
                json = response.json()['rcsb_polymer_entity_feature'][0]
                values = get_json_features(json, keys)
                keys = ['assignment_' + x for x in keys] # two fields are named 'name' so pre-pended with _assignment
                polymer_entity_dict.update(dict(zip(keys, values)))
            except:
                values = [None for key in keys]
                keys = ['assignment_' + x for x in keys] # two fields are named 'name' so pre-pended with _assignment
                polymer_entity_dict.update(dict(zip(keys, values)))
            
            json = response.json()['rcsb_polymer_entity']['rcsb_macromolecular_names_combined'][0]
            keys = ['name']
            values = get_json_features(json, keys)
            
            keys = ['rcsb_macromolecular_names_combined_' + x for x in keys] # two fields are named 'name' so pre-pended with _assignment
            polymer_entity_dict.update(dict(zip(keys, values)))
            
            keys = ['reference_database_name', 'reference_database_accession']
            try:
                json = response.json()['rcsb_polymer_entity_align'][0]
                values = get_json_features(json, keys)
                polymer_entity_dict.update(dict(zip(keys, values)))
            except:
                values = [None for key in keys]
                polymer_entity_dict.update(dict(zip(keys, values)))

            df = df.append(polymer_entity_dict, ignore_index = True)

    return df

def generate_stem_list(file_paths):
    stem_list = []
    for i in range(len(file_paths)):
        path = file_paths[i]
        stem = path.stem
        if stem not in stem_list:
            stem_list.append(stem)

    return stem_list       

def get_json_features(json, keys):
    """Function to get json features

    Args:
        json: json file
        keys: List of keys of the features wanted

    Yields:
        feature_list: values from keys
    """
    feature_list = []
    for key in keys:
        try:
            feature = json[key]
        except:
            feature = ''
            
        feature_list.append(feature)
            
    return feature_list