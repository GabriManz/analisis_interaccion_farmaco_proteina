from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np

# Crear cliente para moléculas, actividades y objetivos
molecule = new_client.molecule
activity = new_client.activity
target = new_client.target

# Cargar el dataset (cambia 'ruta/del/archivo' a la ubicación de tu archivo)
dataset = pd.read_csv('/content/drive/MyDrive/Bioinformática/TFM/Filtered_Extract_Ligands_Uniprot.csv')

# Añadir una columna 'pChEMBL_specific' para los valores específicos del complejo
dataset['pChEMBL'] = None

# Iterar sobre las filas del dataframe
for index, row in dataset.iterrows():
    inchi_key = row['inchi_key']
    uniprot_id = row['UNIPROT_ID']

    # Verificar si tanto el InChI Key como el UniProt ID no son NaN
    if pd.notna(inchi_key) and pd.notna(uniprot_id):
        try:
            # Buscar el ChEMBL ID del ligando a partir del InChI Key
            mol = molecule.filter(molecule_structures__standard_inchi_key=inchi_key).only(['molecule_chembl_id'])
            if mol:
                chembl_id = mol[0]['molecule_chembl_id']

                # Buscar el ChEMBL ID del objetivo usando el UniProt ID
                target_res = target.filter(target_components__accession=uniprot_id).only(['target_chembl_id'])
                if target_res:
                    target_chembl_id = target_res[0]['target_chembl_id']

                    # Buscar actividades específicas para el par (molecule_chembl_id, target_chembl_id) con valores de pChEMBL
                    res = activity.filter(molecule_chembl_id=chembl_id, target_chembl_id=target_chembl_id, pchembl_value__isnull=False)

                    # Obtener los valores pChEMBL
                    pchembl_values = [float(record['pchembl_value']) for record in res if record['pchembl_value'] is not None]

                    if pchembl_values:
                        # Calcular la mediana de los valores pChEMBL específicos del complejo
                        pchembl_median = np.median(pchembl_values)
                        dataset.at[index, 'pChEMBL'] = pchembl_median
                    else:
                        # Si no hay valores pChEMBL específicos, asignar NaN
                        dataset.at[index, 'pChEMBL'] = np.nan
                else:
                    # Si no se encuentra target_chembl_id para el UniProt ID, asignar NaN
                    dataset.at[index, 'pChEMBL'] = np.nan
            else:
                # Si no se encuentra un ChEMBL ID para el InChI Key, asignar NaN
                dataset.at[index, 'pChEMBL'] = np.nan
        except Exception as e:
            # En caso de error, asignar NaN y mostrar el error
            print(f"Error en la fila {index} con InChI Key: {inchi_key} y UniProt ID: {uniprot_id}, Error: {e}")
            dataset.at[index, 'pChEMBL'] = np.nan

# Cargar Dataset de PDB BIND
pdbbind_data = pd.read_excel('./data/input/pdb_bind_data.xlsx')

# Escoger las columnas que nos interesan del dataset
pdbbind_data = pdbbind_data[["PDB code", "Ligand Name", "pKd pKi pIC50"]]

# Renombrar Columnas
pdbbind_data = pdbbind_data.rename(columns={"PDB code": "Complejos PDB_ID", "Ligand Name": "Ligand PDB_ID", "pKd pKi pIC50": "pChEMBL"})

# Realizar el merge usando 'Complejos PDB_ID' y 'Ligand PDB_ID' como claves
dataset = dataset.merge(pdbbind_data[['Complejos PDB_ID', 'Ligand PDB_ID', 'pChEMBL']], 
                        how="left", 
                        on=["Complejos PDB_ID", "Ligand PDB_ID"])

# Combinar pChEMBL_x y pChEMBL_y en una sola columna llamada pChEMBL
dataset["pChEMBL"] = dataset["pChEMBL_y"].combine_first(dataset["pChEMBL_x"])

# Eliminar las columnas pChEMBL_x y pChEMBL_y
dataset = dataset.drop(columns=["pChEMBL_x", "pChEMBL_y"])

# Contar valores no nulos en la columna 'pChEMBL_specific'
num_pchembl_calculated = dataset['pChEMBL'].notna().sum()
print(f"Número de valores pChEMBL específicos calculados: {num_pchembl_calculated}")
print(f"Total de proteínas en el dataset: {len(dataset)}")

# Pasar el dataset a CSV
dataset.to_csv('./data/processed/dataset_with_pchembl.csv', index=False)
