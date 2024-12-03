"""
##########################
### Documentación del Dataset Generado
##########################

Descripción General:
Este script procesa datos provenientes de ChEMBL, BindingDB y PDBBind para generar un dataset integrado con valores pChEMBL,
pKi, y otras métricas de afinidad ligando-proteína. El dataset final permite el análisis y modelado de interacciones
químico-biológicas mediante valores estandarizados de actividad.

Origen de los Datos:
1. ChEMBL:
   - Valores de pChEMBL asociados a combinaciones de ligandos (InChI Key) y proteínas (UniProt ID).
   - Obtenido mediante consultas a la API de ChEMBL.
2. BindingDB:
   - Valores de afinidad (Ki) extraídos del archivo `BindingDB_All.tsv`.
   - Convertidos a valores pKi para su estandarización.
3. PDBBind:
   - Valores pChEMBL extraídos del archivo `pdb_bind_data.xlsx`.
   - Incluye identificadores de complejos PDB y ligandos asociados.

Objetivo del Dataset:
- Proporcionar un dataset estandarizado para el análisis de interacciones ligando-proteína.
- Integrar datos de múltiples fuentes para facilitar estudios comparativos.
- Generar valores de actividad (pChEMBL) para ligandos y proteínas con datos faltantes imputados mediante valores aleatorios entre cuartiles (Q1, Q3).

Estructura del Dataset Final:
- Columnas principales:
  1. `InChI`: Identificador químico del ligando.
  2. `UNIPROT_ID`: Identificador UniProt de la proteína objetivo.
  3. `pChEMBL`: Valor de actividad calculado o imputado, basado en valores de pChEMBL, pKi o valores faltantes rellenados.
  4. `Complejos PDB_ID`: Identificador PDB del complejo asociado.

Procesamiento:
1. **Obtención de datos desde ChEMBL:**
   - Se realiza una búsqueda por InChI Key y UniProt ID.
   - Se extraen valores pChEMBL mediante consultas a la API con control de tasa y reintentos.
2. **Procesamiento de BindingDB:**
   - Se calculan valores pKi a partir de Ki (nM) con la fórmula `-log10(Ki/1e9)`.
   - Se agrupan los datos por UniProt ID e InChI Key para calcular la mediana de pKi.
3. **Integración de PDBBind:**
   - Se extraen valores pChEMBL asociados a identificadores PDB.
4. **Combinación de valores:**
   - Se priorizan los valores pChEMBL de PDBBind, seguidos por los de ChEMBL y luego pKi.
   - Valores faltantes de pChEMBL se imputan con valores aleatorios entre Q1 y Q3.
5. **Agrupación final:**
   - Se agrupan los datos por `InChI` y `UNIPROT_ID`, conservando la media de `pChEMBL` y el primer `Complejos PDB_ID`.

Archivos Generados:
1. Dataset Final:
   - Archivo de salida: `dataset_with_pchembl.csv`.
   - Contiene los valores de pChEMBL integrados y procesados para cada combinación única de ligando y proteína.

Ubicaciones de Archivos:
- Archivos de entrada:
  1. `Filtered_Extract_Ligands_Uniprot.csv`: Ligandos y proteínas originales.
  2. `BindingDB_All.tsv`: Datos de afinidad de BindingDB.
  3. `pdb_bind_data.xlsx`: Datos de PDBBind.
- Archivo de salida:
  - `./data/processed/dataset_with_pchembl.csv`.

Limitaciones:
- La calidad de los valores pChEMBL depende de la disponibilidad y precisión de los datos en ChEMBL, BindingDB y PDBBind.
- Ligandos y proteínas sin valores conocidos tendrán imputaciones aleatorias dentro de un rango estadístico.

Tamaño Aproximado:
- Dependerá del número de combinaciones de ligandos y proteínas disponibles en los archivos de entrada.
"""


from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
import time
from random import uniform

# Crear cliente para ChEMBL
molecule = new_client.molecule
activity = new_client.activity
target = new_client.target

# Ruta de los archivos
original_dataset_path = './data/processed/Filtered_Extract_Ligands_Uniprot.csv'
bindingdb_file_path = './data/input/BindingDB_All.tsv'
pdbbind_file_path = './data/input/pdb_bind_data.xlsx'
output_dataset_path = './data/processed/dataset_with_pchembl.csv'

# Cargar el dataset principal
original_dataset = pd.read_csv(original_dataset_path)
original_dataset['pChEMBL'] = None

# Función para obtener datos de ChEMBL con reintentos y control de tasa
def fetch_chembl_data_with_rate_limiting(inchi_key, uniprot_id, retries=3, base_delay=5):
    """Función para obtener datos de ChEMBL con reintentos y control de tasa."""
    for attempt in range(1, retries + 1):
        try:
            # Buscar ChEMBL ID del ligando
            mol = molecule.filter(molecule_structures__standard_inchi_key=inchi_key).only(['molecule_chembl_id'])
            mol_list = list(mol)
            if mol_list:
                chembl_id = mol_list[0]['molecule_chembl_id']

                # Buscar ChEMBL ID del target
                target_res = target.filter(target_components__accession=uniprot_id).only(['target_chembl_id'])
                target_list = list(target_res)
                if target_list:
                    target_chembl_id = target_list[0]['target_chembl_id']

                    # Buscar actividades específicas con valores de pChEMBL
                    res = activity.filter(
                        molecule_chembl_id=chembl_id,
                        target_chembl_id=target_chembl_id,
                        pchembl_value__isnull=False
                    )
                    res_list = list(res)
                    if res_list:
                        pchembl_values = [
                            float(record['pchembl_value']) for record in res_list if record.get('pchembl_value') is not None
                        ]
                        if pchembl_values:
                            return np.median(pchembl_values)  # Retorna la mediana
            return None  # Si no se encuentra nada
        except Exception as e:
            print(f"Error para {inchi_key}, {uniprot_id} en intento {attempt}/{retries}: {e}")
            time.sleep(base_delay * attempt + uniform(1, 3))  # Espera exponencial con variabilidad
    return None

# Buscar combinaciones únicas de UNIPROT_ID e inchi_key
unique_combinations = original_dataset[['inchi_key', 'UNIPROT_ID']].drop_duplicates()

# Diccionario para almacenar resultados de ChEMBL
pchembl_dict = {}

# Procesar cada combinación única
for index, row in unique_combinations.iterrows():
    inchi_key = row['inchi_key']
    uniprot_id = row['UNIPROT_ID']
    if pd.notna(inchi_key) and pd.notna(uniprot_id):
        pchembl_value = fetch_chembl_data_with_rate_limiting(inchi_key, uniprot_id)
        if pchembl_value is not None:
            pchembl_dict[(inchi_key, uniprot_id)] = pchembl_value

# Asignar valores de pChEMBL obtenidos al dataset original
original_dataset['pChEMBL'] = original_dataset.apply(
    lambda x: pchembl_dict.get((x['inchi_key'], x['UNIPROT_ID']), x['pChEMBL']),
    axis=1
)

# Procesar y agrupar datos de BindingDB
columns_to_load = [
    'Ligand InChI Key',
    'UniProt (SwissProt) Primary ID of Target Chain',
    'Ki (nM)'
]

# Cargar BindingDB en fragmentos
chunksize = 100000
bindingdb_chunks = []
for chunk in pd.read_csv(bindingdb_file_path, sep="\t", usecols=columns_to_load, chunksize=chunksize, on_bad_lines="skip", low_memory=True):
    bindingdb_chunks.append(chunk)

# Combinar todos los fragmentos en un DataFrame final
bindingdb_data = pd.concat(bindingdb_chunks, ignore_index=True)

# Limpiar y calcular pKi
bindingdb_data['Ki (nM)'] = pd.to_numeric(bindingdb_data['Ki (nM)'], errors='coerce')
bindingdb_data.loc[bindingdb_data['Ki (nM)'] == 0, 'Ki (nM)'] = np.nan
bindingdb_data['pKi'] = bindingdb_data['Ki (nM)'].apply(
    lambda x: -np.log10(x / 1e9) if pd.notna(x) else np.nan
)

# Agrupar por UniProt ID e InChI Key y calcular la mediana de pKi
bindingdb_grouped = bindingdb_data.groupby(
    ['UniProt (SwissProt) Primary ID of Target Chain', 'Ligand InChI Key']
)['pKi'].median().reset_index()
bindingdb_grouped.rename(columns={"pKi": "median_pKi"}, inplace=True)

# Realizar el merge entre el dataset original y BindingDB
merged_dataset = pd.merge(
    original_dataset,
    bindingdb_grouped,
    how="left",
    left_on=["UNIPROT_ID", "inchi_key"],
    right_on=["UniProt (SwissProt) Primary ID of Target Chain", "Ligand InChI Key"]
)

# Integrar datos de PDBBind
pdbbind_data = pd.read_excel(pdbbind_file_path)
pdbbind_data = pdbbind_data.rename(columns={
    "PDB code": "Complejos PDB_ID",
    "Ligand Name": "Ligand PDB_ID",
    "pKd pKi pIC50": "pChEMBL"
})
merged_dataset = merged_dataset.merge(
    pdbbind_data[['Complejos PDB_ID', 'Ligand PDB_ID', 'pChEMBL']],
    how="left",
    on=["Complejos PDB_ID", "Ligand PDB_ID"]
)

# Combinar los valores de pChEMBL y pKi priorizando pChEMBL
merged_dataset["final_pChEMBL"] = merged_dataset["pChEMBL_y"].combine_first(
    merged_dataset["pChEMBL_x"].combine_first(merged_dataset["median_pKi"])
)

# Eliminar columnas auxiliares
merged_dataset.drop(columns=["pChEMBL_x", "pChEMBL_y", "median_pKi", "UniProt (SwissProt) Primary ID of Target Chain", "Ligand InChI Key"], inplace=True)

# Agrupar por 'InChI' y 'UNIPROT_ID' para calcular la media de 'final_pChEMBL' y conservar el primer 'Complejos PDB_ID'
mean_values = (
    merged_dataset
    .groupby(['InChI', 'UNIPROT_ID'], as_index=False)
    .agg({
        'final_pChEMBL': 'mean',           # Calcular la media de 'final_pChEMBL'
        'Complejos PDB_ID': 'first'       # Conservar el primer valor de 'Complejos PDB_ID'
    })
    .rename(columns={'final_pChEMBL': 'pChEMBL'})  # Renombrar la columna
)

# Calcular los cuartiles (Q1 y Q3) de 'pChEMBL'
q1, q3 = mean_values['pChEMBL'].quantile([0.25, 0.75])

# Rellenar valores NaN en 'pChEMBL' con valores aleatorios entre Q1 y Q3
mean_values['pChEMBL'] = mean_values['pChEMBL'].apply(
    lambda x: np.random.uniform(q1, q3) if pd.isna(x) else x
)

# Guardar el dataset final
mean_values.to_csv(output_dataset_path, index=False)

print(f"Dataset procesado y guardado en: {output_dataset_path}")