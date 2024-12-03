"""
##########################
### Documentación del Dataset
##########################

Descripción General:
Este dataset combina información de ligandos y decoys procesados para tareas de modelado químico-biológico.
Se incluye información química, representaciones estructurales y métricas de similitud, junto con datos funcionales
como representaciones n-hot de los grupos funcionales más frecuentes.

Origen de los Datos:
1. Ligandos iniciales:
   - Archivo de entrada: `filtered_extract_ligands_uniprot.csv`.
   - Contiene identificadores de ligandos en formato InChI y metadatos asociados.
   - Convertidos a formato SMILES para su procesamiento.
2. Decoys:
   - Carpeta de entrada: `DUD-E_results/`.
   - Archivos `.picked` con información sobre ligandos y decoys asociados.
3. Representaciones funcionales:
   - Archivo intermedio: `fgs_decoy.csv`.
   - Incluye pseudo-SMILES generados a partir de las estructuras de los decoys.

Objetivo del Dataset:
- Facilitar estudios comparativos entre ligandos activos y sus decoys (compuestos señuelo).
- Permitir análisis de similitud estructural mediante la métrica Tanimoto.
- Proporcionar una base para modelar interacciones químico-biológicas mediante representaciones funcionales.

Estructura del Dataset:
1. Ligandos y Decoys:
   - `Ligand_SMILES`: Representación SMILES del ligando original.
   - `Ligand_ID`: Identificador único del ligando.
   - `Decoy_SMILES`: Representación SMILES del decoy asociado.
   - `Tanimoto_Similarity`: Similitud Tanimoto entre el ligando y el decoy.
   - `Actividad`: Valor binario (0 para decoys, 1 para ligandos activos).
2. Representaciones funcionales (n-hot):
   - `n_hot_fgs`: Representación n-hot de los grupos funcionales más frecuentes.

Archivos Generados:
1. `valid_smiles.txt`: Lista de SMILES válidos extraídos de los ligandos iniciales.
2. `ligand_decoy.csv`: Dataset combinado de ligandos y decoys, incluyendo similitudes Tanimoto.
3. `fgs_decoy.csv`: Representación funcional en formato pseudo-SMILES.
4. `fgs_decoy_n_hot.csv`: Representaciones n-hot de los grupos funcionales más frecuentes.

Ubicaciones de Archivos:
- Todos los archivos generados se guardan en el directorio `./data/processed/`.

Tamaño Aproximado:
- Depende del número de ligandos y decoys procesados. Los datasets suelen contener:
  - Miles de ligandos iniciales.
  - Decoys generados para cada ligando (dependiendo del algoritmo DUD-E).
  - Grupos funcionales representados en las columnas n-hot.
"""


import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from fgs import generate_fgs_decoy

# 1. Leer el archivo CSV inicial
def load_data(file_path):
    return pd.read_csv(file_path)

# 2. Convertir InChI a SMILES
def inchi_to_smiles(inchi):
    try:
        mol = Chem.MolFromInchi(inchi)  # Convertir InChI a molécula
        if mol:
            return Chem.MolToSmiles(mol)  # Convertir la molécula a SMILES
    except Exception:
        return None

# 3. Filtrar SMILES válidos
def filter_valid_smiles(data):
    data['SMILES'] = data['InChI'].apply(inchi_to_smiles)
    return data['SMILES'].drop_duplicates().dropna()

# 4. Calcular la similitud de Tanimoto
def tanimoto_similarity(smiles1, smiles2):
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 and mol2:
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            return DataStructs.TanimotoSimilarity(fp1, fp2)
    except Exception:
        return None

# 5. Procesar los archivos `.picked` de una carpeta
def process_decoy_folder(folder_path):
    ligands_and_decoys = []
    for file in os.listdir(folder_path):
        if file.endswith(".picked"):
            with open(os.path.join(folder_path, file), 'r') as f:
                lines = f.readlines()
                if not lines:
                    continue
                # Primera línea: Ligando original
                ligand_info = lines[0].strip().split('\t')
                ligand_smiles = ligand_info[1]
                ligand_id = ligand_info[3]

                # Resto de líneas: Decoys
                for line in lines[1:]:
                    decoy_info = line.strip().split('\t')
                    decoy_smiles = decoy_info[0]
                    ligands_and_decoys.append({
                        'Ligand_SMILES': ligand_smiles,
                        'Ligand_ID': ligand_id,
                        'Decoy_SMILES': decoy_smiles
                    })
    return pd.DataFrame(ligands_and_decoys)

# 6. Procesar todas las carpetas de decoys
def process_all_decoy_folders(parent_folder):
    all_decoys = []
    for chunk_folder in os.listdir(parent_folder):
        chunk_path = os.path.join(parent_folder, chunk_folder)
        if os.path.isdir(chunk_path):
            print(f"Procesando carpeta: {chunk_folder}")
            chunk_data = process_decoy_folder(chunk_path)
            all_decoys.append(chunk_data)
    return pd.concat(all_decoys, ignore_index=True)

# 7. Filtrar los decoys menos similares
def find_least_similar_decoys(data, top_n=20):
    filtered_data = []
    for ligand_id, ligand_data in data.groupby('Ligand_ID'):
        ligand_smiles = ligand_data['Ligand_SMILES'].iloc[0]
        ligand_data['Tanimoto_Similarity'] = ligand_data['Decoy_SMILES'].apply(
            lambda decoy: tanimoto_similarity(ligand_smiles, decoy)
        )
        least_similar_decoys = ligand_data.nsmallest(top_n, 'Tanimoto_Similarity')
        filtered_data.append(least_similar_decoys)
    return pd.concat(filtered_data, ignore_index=True)

# 9. Generar representaciones "n-hot" de los grupos funcionales
def generate_n_hot_representation(input_path, output_path):
    """
    Genera una representación n-hot de los grupos funcionales más frecuentes.
    
    Parámetros:
    - input_path: Ruta del archivo CSV de entrada con pseudo-SMILES.
    - output_path: Ruta del archivo CSV de salida con las representaciones n-hot.
    """
    print(f"Procesando {input_path} para generar representaciones n-hot...")
    
    # Cargar los datos desde el archivo
    data = pd.read_csv(input_path)
    
    # Contar la frecuencia de cada grupo funcional en la columna 'psmis'
    group_counts = data['psmis'].value_counts()
    
    # Seleccionar los 100 grupos funcionales más frecuentes
    top_100_groups = group_counts.head(100).index
    
    # Crear una columna con listas de grupos funcionales
    data['functional_groups'] = data['psmis'].apply(lambda x: [x])
    
    # Agrupar por 'Decoy_SMILES' y consolidar los grupos funcionales en una lista única
    grouped_data = data.groupby('Decoy_SMILES').agg({
    'functional_groups': 'sum',  
    'Ligand_SMILES': 'first',   
    'Actividad': 'first'         
    }).reset_index()
    
    # Contar la ocurrencia de cada grupo funcional en las listas
    grouped_data['n_hot_fgs'] = grouped_data['functional_groups'].apply(
        lambda x: [x.count(group) for group in top_100_groups]
    )
    
    # Guardar los resultados en un archivo CSV
    grouped_data[['Ligand_SMILES', 'Decoy_SMILES', 'Actividad', 'n_hot_fgs']].to_csv(output_path, index=False)
    print(f"Representaciones n-hot guardadas en: {output_path}")

# 10. Pipeline principal actualizado
def main():
    # Archivos y carpetas
    inchi_file_path = './data/processed/filtered_extract_ligands_uniprot.csv'
    parent_folder_path = './data/DUD-E_results/'
    output_file_path = './data/processed/ligand_decoy.csv'
    pseudo_smiles_output_path = './data/processed/fgs_decoy.csv'
    n_hot_output_path = './data/processed/fgs_decoy_n_hot.csv'

    # Paso 1: Leer archivo y convertir InChI a SMILES
    print("Cargando datos y convirtiendo InChI a SMILES...")
    data = load_data(inchi_file_path)
    data_smiles = filter_valid_smiles(data)

    # Guardar los SMILES válidos para verificar
    data_smiles.to_csv('./data/processed/valid_smiles.txt', index=False, header=False)

    # Paso 2: Procesar carpetas de decoys
    print("Procesando carpetas de decoys...")
    decoy_data = process_all_decoy_folders(parent_folder_path)

    # Paso 3: Filtrar decoys menos similares
    print("Filtrando los decoys menos similares...")
    filtered_decoy_data = find_least_similar_decoys(decoy_data, top_n=20)

    # Paso 4: Añadir columna de actividad
    filtered_decoy_data["Actividad"] = 0

    # Paso 5: Guardar resultados intermedios
    print(f"Guardando resultados intermedios en {output_file_path}...")
    filtered_decoy_data.to_csv(output_file_path, index=False)

    # Paso 6: Generar pseudo-SMILES
    print("Generando pseudo-SMILES...")
    generate_fgs_decoy(output_file_path, pseudo_smiles_output_path)

    # Paso 7: Generar representaciones n-hot
    print("Generando representaciones n-hot de grupos funcionales...")
    generate_n_hot_representation(pseudo_smiles_output_path, n_hot_output_path)

    print("Proceso completado.")

# Ejecutar el pipeline
if __name__ == "__main__":
    main()
