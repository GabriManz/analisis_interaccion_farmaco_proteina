# -*- coding: utf-8 -*-
"""
Script: extract_ligands_uniprot.py

Descripción:
Este script permite identificar ligandos asociados a estructuras de proteínas descargadas desde el PDB y
asociarlas a identificadores de UniProt. 

Flujo principal:
1. Obtener el UniProt ID de cada PDB ID utilizando la API de UniProt.
2. Extraer y agregar ligandos a partir de los archivos CIF de proteínas descargadas.
3. Generar un archivo que asocia los complejos PDB y ligandos con sus UniProt IDs.

Salidas:
- CSV final con asociaciones de complejos PDB, ligandos, y UniProt IDs.
"""

import requests
import time
import pandas as pd
import os
import glob
from Bio.PDB.MMCIFParser import MMCIFParser

# Función para obtener UniProt ID desde un PDB ID usando la API de UniProt
def get_uniprot_id_from_pdb(pdb_id):
    url = 'https://rest.uniprot.org/idmapping/run'
    files = {
        'ids': (None, pdb_id),
        'from': (None, 'PDB'),
        'to': (None, 'UniProtKB'),
    }

    try:
        response = requests.post(url, files=files)
        response.raise_for_status()
        job_id = response.json().get('jobId')
    except requests.RequestException as e:
        print(f'Error en la solicitud: {e}')
        return [(pdb_id, None)]

    status_url = f'https://rest.uniprot.org/idmapping/status/{job_id}'

    while True:
        try:
            response = requests.get(status_url)
            response.raise_for_status()
            status_info = response.json()
            if 'results' in status_info:
                result_data = [(result['from'], result['to']['primaryAccession']) for result in status_info['results']]
                return result_data if result_data else [(pdb_id, None)]
            else:
                print('El trabajo aún está en progreso. Esperando...')
                time.sleep(10)
        except requests.RequestException as e:
            print(f'Error en la solicitud: {e}')
            return [(pdb_id, None)]

# Función para extraer y agregar ligandos a partir de archivos CIF
def extract_and_add_ligands(df_merged, cif_directory):
    parser = MMCIFParser(QUIET=True)
    ligands_dict = {}

    pdb_ids = df_merged['Complejos PDB_ID'].unique()
    for pdb_id in pdb_ids:
        cif_file = os.path.join(cif_directory, f'{pdb_id}.cif')
        try:
            structure = parser.get_structure(pdb_id, cif_file)
            hetatm_list = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != ' ':
                            hetatm_info = f"{residue.resname} {chain.id}"
                            if hetatm_info not in hetatm_list:
                                hetatm_list.append(hetatm_info)
            ligands_dict[pdb_id] = ', '.join(hetatm_list)
        except FileNotFoundError:
            print(f"Archivo CIF no encontrado para PDB_ID: {pdb_id}")
            ligands_dict[pdb_id] = None

    df_merged['All Ligands'] = df_merged['Complejos PDB_ID'].map(ligands_dict)
    return df_merged

# Función para procesar un lote de archivos CIF y obtener UniProt IDs
def process_batch(batch):
    results = []
    for archivo in batch:
        pdb_id = os.path.basename(archivo).split('.')[0]
        data = get_uniprot_id_from_pdb(pdb_id)
        if data:
            results.extend(data)
    return results

# Función principal que ejecuta el flujo de trabajo completo
def main_extract_ligands_uniprot():
    # Rutas relativas a la estructura del proyecto
    cif_directory = "./data/pdb_files/"
    resultados_path = "./data/processed/result_pdb.csv"
    blacklist_path = "./data/input/blacklist.txt"
    output_path = "./data/processed/Filtered_Extract_Ligands_Uniprot.csv"

    # Procesar archivos CIF en lotes y obtener UniProt IDs
    archivos_cif = glob.glob(f"{cif_directory}/*.cif")
    results_list = []
    batch_size = 500

    for i in range(0, len(archivos_cif), batch_size):
        batch = archivos_cif[i:i + batch_size]
        results = process_batch(batch)
        results_list.extend(results)
        print(f'Processed batch {i // batch_size + 1} of {len(archivos_cif) // batch_size + 1}')

    # Crear DataFrame con los resultados y limpiar los datos
    df_uniprot = pd.DataFrame(results_list, columns=['Complejos PDB_ID', 'UNIPROT_ID'])
    df_uniprot['Complejos PDB_ID'] = df_uniprot['Complejos PDB_ID'].astype(str)
    df_uniprot = df_uniprot.groupby('Complejos PDB_ID')['UNIPROT_ID'].apply(lambda x: ','.join(filter(None, x))).reset_index()

    # Leer resultados anteriores y expandir
    resultados = pd.read_csv(resultados_path, dtype={'Complejos PDB_ID': str})
    resultados_expanded = resultados.assign(Complejos_PDB_ID=resultados['Complejos PDB_ID'].str.split()).explode('Complejos_PDB_ID')
    resultados_expanded.rename(columns={'Complejos_PDB_ID': 'Complejos PDB_ID'}, inplace=True)

    # Merge de datos de UniProt con datos anteriores
    df_merged = df_uniprot.merge(resultados_expanded[['Complejos PDB_ID', 'Ligand PDB_ID', 'InChI', 'inchi_key']], on='Complejos PDB_ID', how='left')

    # Agregar ligandos desde archivos CIF
    df_merged = extract_and_add_ligands(df_merged, cif_directory)

    # Filtrar usando blacklist
    with open(blacklist_path, 'r') as file:
        black_list = [line.split(', ')[0].strip() for line in file if ', ' in line]

    num_ligands_before = len(df_merged)
    df_merged_filtered = df_merged[~df_merged['Ligand PDB_ID'].isin(black_list)]
    num_ligands_after = len(df_merged_filtered)

    print(f'Número de ligandos eliminados por blacklist: {num_ligands_before - num_ligands_after}')

    # Guardar DataFrame final en CSV
    df_merged_filtered.to_csv(output_path, index=False)
    print(f"Archivo guardado en: {output_path}")
