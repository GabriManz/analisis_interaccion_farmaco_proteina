# -*- coding: utf-8 -*-

"""
Script: interacciones_proteina.py

Descripción:
Este script procesa un conjunto de fármacos orales obtenidos de DrugBank, filtrando y estandarizando las moléculas.
A continuación, mapea los fármacos a identificadores PDB para obtener sus complejos proteína-ligando asociados,
y guarda los resultados en archivos para su uso posterior.

Pasos principales:
1. Carga y limpieza de datos.
2. Estandarización y filtrado de estructuras moleculares.
3. Búsqueda de identificadores PDB para los fármacos en el Protein Data Bank.
4. Obtención de complejos proteína-ligando asociados a los fármacos.
5. Exportación de los resultados a archivos CSV y TXT.
"""

# Librerías
import pandas as pd
from rdkit import Chem
from chembl_structure_pipeline import standardizer as sdz

# Cargar datos de fármacos orales
def cargar_datos(file_path):
    """Cargar el archivo CSV en un DataFrame y obtener el tamaño inicial para referencia."""
    df = pd.read_csv(file_path, delimiter=";")
    tamaño_inicial = df.shape[0]
    print(f"Tamaño inicial del dataset: {tamaño_inicial} compuestos.")
    return df

# Estandarización y filtrado de datos
def estandarizar_y_filtrar(df):
    """Estandariza y filtra las moléculas, eliminando duplicados y compuestos no válidos."""
    # Generar las moléculas a partir de la cadena InChI
    df["mol"] = df.inchi.apply(Chem.MolFromInchi)

    # Obtener la molécula "padre" estandarizada para cada molécula
    df["pmol"] = df.mol.apply(lambda x: sdz.get_parent_mol(sdz.standardize_mol(x))[0])

    # Filtrar compuestos con varias moléculas (excluyendo aquellos con múltiples componentes en SMILES)
    df = df[df.pmol.apply(lambda x: "." not in Chem.MolToSmiles(x))]

    # Eliminar compuestos duplicados utilizando el InChIKey de la molécula estandarizada
    df["inchi_key"] = df.pmol.apply(lambda x: Chem.MolToInchiKey(x))
    df = df.drop_duplicates(subset=["inchi_key"])

    # Filtrar los compuestos que tienen más de 6 carbonos
    df = df[df.pmol.apply(lambda x: sum(1 for atom in x.GetAtoms() if atom.GetSymbol() == 'C') > 6)]
    
    print(f"Después del filtrado, quedan {df.shape[0]} compuestos.")
    return df

# Guardar resultados preliminares de DrugBank
def guardar_drugbank(df, output_path):
    """Guarda un DataFrame con información seleccionada de los fármacos en un archivo CSV."""
    result_drugbank = df[['hmdb_id', 'name', 'inchi_key', 'inchi']]
    result_drugbank.to_csv(output_path, index=False)
    print(f"Archivo de DrugBank guardado en: {output_path}")

# Búsqueda de identificadores PDB
def buscar_pdb_ids(df):
    """Busca PDB IDs correspondientes a los fármacos utilizando InChI y InChIKey."""
    # Leer el archivo de InChI
    inchi_file_url = "http://ligand-expo.rcsb.org/dictionaries/Components-inchi.ich"
    inchi_data = pd.read_csv(inchi_file_url, sep='\t', header=None)
    inchi_data.columns = ['InChI', 'Ligand PDB_ID', 'Name']

    # Filtrar entradas que coinciden con los InChIs en el dataset
    matching_inchis = inchi_data[inchi_data['InChI'].isin(set(df.inchi))]

    # Leer el archivo de InChI-Key
    inchi_key_url = "http://ligand-expo.rcsb.org/dictionaries/Components-inchikey.ich"
    inchi_key_data = pd.read_csv(inchi_key_url, sep='\t', header=None)
    inchi_key_data.columns = ['inchi_key', 'Ligand PDB_ID', 'Name']

    # Filtrar entradas que coinciden con los InChI-Keys en el dataset
    matching_inchis_keys = inchi_key_data[inchi_key_data['inchi_key'].isin(set(df['inchi_key']))]

    # Unir los resultados de InChI e InChI-Key
    result = pd.merge(matching_inchis, matching_inchis_keys, on=['Ligand PDB_ID', 'Name'])
    print(f"Se encontraron {result.shape[0]} compuestos con PDB IDs correspondientes.")
    return result

# Búsqueda de complejos proteína-ligando
def buscar_complejos_ligando(result):
    """Busca los complejos proteína-ligando asociados a los PDB IDs encontrados."""
    complejo_ligando_url = "http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd"
    complejo_ligando = pd.read_csv(complejo_ligando_url, sep='\t', header=None)
    complejo_ligando.columns = ['Ligand PDB_ID', 'Complejos PDB_ID']

    # Realizar el merge para obtener los complejos proteína-ligando
    result_pdb = pd.merge(result, complejo_ligando, how='left', on='Ligand PDB_ID')
    return result_pdb

# Guardar resultados finales
def guardar_resultados(result_pdb, csv_path, txt_path):
    """Guarda el DataFrame con los complejos en CSV y exporta una lista única de PDB IDs en TXT."""
    result_pdb.to_csv(csv_path, index=False)
    print(f"Archivo CSV de resultados guardado en: {csv_path}")

    # Extraer y guardar IDs únicos de los complejos en un archivo de texto
    pdb_ids_unicos = set()
    for ids in result_pdb['Complejos PDB_ID']:
        pdb_ids_unicos.update(str(ids).split())

    with open(txt_path, 'w') as file:
        for pdb_id in pdb_ids_unicos:
            file.write(pdb_id + '\n')
    print(f"Archivo TXT de PDB IDs únicos guardado en: {txt_path}")

# Función principal que coordina todo el proceso
def main_interacciones_proteina():
    # Definir las rutas de entrada y salida
    file_path = "./data/input/drugs.csv"
    drugbank_output_path = "./data/processed/result_drugbank.csv"
    result_csv_path = "./data/processed/result_pdb.csv"
    result_txt_path = "./data/processed/complejos_PDB.txt"

    print("Cargando y estandarizando datos de fármacos...")
    df = cargar_datos(file_path)
    df = estandarizar_y_filtrar(df)

    print("Guardando resultados preliminares de DrugBank...")
    guardar_drugbank(df, drugbank_output_path)

    print("Buscando identificadores PDB...")
    result = buscar_pdb_ids(df)

    print("Buscando complejos proteína-ligando...")
    result_pdb = buscar_complejos_ligando(result)

    print("Guardando resultados finales...")
    guardar_resultados(result_pdb, result_csv_path, result_txt_path)
    print("Proceso completado.")
