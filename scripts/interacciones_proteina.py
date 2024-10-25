# -*- coding: utf-8 -*-
"""
Script: interacciones_proteina.py
@author: Gabriel

Descripción:
Este script procesa un archivo de datos de fármacos para generar un dataset principal con moléculas estandarizadas.
El flujo de trabajo incluye los siguientes pasos:
1. Filtrado y estandarización de datos de moléculas.
2. Identificación de moléculas únicas.
3. Búsqueda de identificadores PDB y complejos proteína-ligando.

Salidas:
- Un archivo CSV con el dataset procesado y otro con las interacciones proteína-ligando encontradas en el PDB.
"""

import pandas as pd
from rdkit import Chem
from chembl_structure_pipeline import standardizer as sdz

# Cargar el archivo de datos de fármacos
file_path = "E:/Máster Bioinformática/TFM/drugs.csv"
df = pd.read_csv(file_path, delimiter=";")

# Mostrar el tamaño inicial del dataset
tamaño_inicial = df.shape[0]
print(f"Tamaño inicial del dataset: {tamaño_inicial} registros")

# Generar moléculas RDKit a partir del InChI y estandarizarlas
df["mol"] = df.inchi.apply(Chem.MolFromInchi)
df["pmol"] = df.mol.apply(lambda x: sdz.get_parent_mol(sdz.standardize_mol(x))[0])

# Filtrar compuestos que contengan una sola molécula (sin el caracter ".")
df = df[df.pmol.apply(lambda x: "." not in Chem.MolToSmiles(x))]

# Eliminar compuestos duplicados usando el InChIKey
df["inchi_key"] = df.pmol.apply(lambda x: Chem.MolToInchiKey(x))
df = df.drop_duplicates(subset=["inchi_key"])
tamaño_despues_duplicados = df.shape[0]
print(f"Compuestos duplicados eliminados: {tamaño_inicial - tamaño_despues_duplicados}")

# Función para filtrar compuestos con más de 6 carbonos
def tiene_mas_de_seis_carbonos(mol):
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    return num_carbons > 6

df = df[df.pmol.apply(tiene_mas_de_seis_carbonos)]
tamaño_despues_filtro = df.shape[0]
print(f"Compuestos con menos de 6 carbonos eliminados: {tamaño_despues_duplicados - tamaño_despues_filtro}")

print(f"Tamaño final del dataset: {tamaño_despues_filtro} registros")

# Guardar el dataset procesado
result_drugbank = df[['hmdb_id', 'name', 'inchi_key', 'inchi']]
result_drugbank.to_csv('result_drugbank.csv', index=False)

# Búsqueda de fármacos en el PDB usando InChI
# Leer datos de InChI y unir con el dataset procesado
inchi_file_url = "http://ligand-expo.rcsb.org/dictionaries/Components-inchi.ich"
inchi_data = pd.read_csv(inchi_file_url, sep='\t', header=None, names=['InChI', 'Ligand PDB_ID', 'Name'])
inchi_set = set(df.inchi)

matching_inchis = inchi_data[inchi_data['InChI'].isin(inchi_set)]
print(f"Fármacos coincidentes en el PDB usando InChI: {matching_inchis.shape[0]}")

# Búsqueda de fármacos en el PDB usando InChIKey
inchi_key_url = "http://ligand-expo.rcsb.org/dictionaries/Components-inchikey.ich"
inchi_key_data = pd.read_csv(inchi_key_url, sep='\t', header=None, names=['inchi_key', 'Ligand PDB_ID', 'Name'])
inchi_key_set = set(df['inchi_key'])

matching_inchi_keys = inchi_key_data[inchi_key_data['inchi_key'].isin(inchi_key_set)]
print(f"Fármacos coincidentes en el PDB usando InChIKey: {matching_inchi_keys.shape[0]}")

# Unir coincidencias por InChI e InChIKey en un solo DataFrame
result = pd.merge(
    matching_inchis,
    matching_inchi_keys,
    on=['Ligand PDB_ID', 'Name'],
    how='outer'
)

# Búsqueda de complejos proteína-ligando en el PDB
complejo_ligando_url = "http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd"
complejo_ligando = pd.read_csv(complejo_ligando_url, sep='\t', header=None, names=['Ligand PDB_ID', 'Complejos PDB_ID'])
result_pdb = pd.merge(result, complejo_ligando, on='Ligand PDB_ID', how='left')

# Guardar el DataFrame de interacciones proteína-ligando
result_pdb.to_csv('result_pdb.csv', index=False)
print("Interacciones proteína-ligando guardadas en 'result_pdb.csv'")

# Crear archivo de texto con IDs de los complejos PDB
pdb_ids = result_pdb['Complejos PDB_ID'].dropna().unique()
file_path_ids = 'complejos_PDB.txt'

with open(file_path_ids, 'w') as file:
    for pdb_id in pdb_ids:
        file.write(f"{pdb_id}\n")

print(f"Archivo de IDs de complejos PDB creado en '{file_path_ids}'")
