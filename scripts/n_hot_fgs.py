"""
##########################
### Documentación del Dataset Generado
##########################

Descripción General:
Este script procesa datos relacionados con grupos funcionales en complejos PDB para generar una representación n-hot. 
La representación n-hot captura la frecuencia de los 100 grupos funcionales más comunes en cada complejo PDB, permitiendo
un análisis estructural y funcional simplificado.

Origen de los Datos:
- Archivo de entrada: `fgs_pdb.csv`.
  - Contiene información sobre los grupos funcionales identificados en los complejos PDB.
  - La columna `pseudo_smiles` incluye representaciones simplificadas de los grupos funcionales.

Objetivo del Dataset:
- Convertir la información de grupos funcionales en complejos PDB a una representación n-hot.
- Facilitar el análisis funcional y estructural basado en grupos funcionales más frecuentes.
- Proporcionar una base para modelado predictivo y análisis de características funcionales.

Estructura del Dataset Final:
- Columnas principales:
  1. `Complejos PDB_ID`: Identificador único del complejo PDB.
  2. `n_hot_fgs`: Vector n-hot que representa la frecuencia de los 100 grupos funcionales más frecuentes en el complejo.

Procesamiento:
1. **Frecuencia de grupos funcionales:**
   - Se calcula la frecuencia de aparición de cada grupo funcional en la columna `pseudo_smiles`.
   - Se seleccionan los 100 grupos funcionales más frecuentes para construir el vector n-hot.
2. **Agrupación por Complejos PDB:**
   - Se agrupan los datos por `Complejos PDB_ID`, consolidando las listas de grupos funcionales para cada complejo.
3. **Cálculo del vector n-hot:**
   - Se genera un vector n-hot para cada complejo, indicando la frecuencia de cada uno de los 100 grupos funcionales seleccionados.
4. **Guardado del dataset:**
   - El resultado final se guarda en un archivo CSV que incluye los identificadores de los complejos PDB y sus vectores n-hot.

Archivos Generados:
1. Dataset Final:
   - Archivo de salida: `fgs_n_hot.csv`.
   - Contiene los vectores n-hot generados para cada complejo PDB.

Ubicaciones de Archivos:
- Archivo de entrada: `./data/processed/fgs_pdb.csv`.
- Archivo de salida: `./data/processed/fgs_n_hot.csv`.

Limitaciones:
- Solo se procesan los 100 grupos funcionales más frecuentes, lo que puede no capturar toda la diversidad funcional.
- La calidad de la representación n-hot depende de la precisión de los datos en `pseudo_smiles`.

Tamaño Aproximado:
- El tamaño del dataset final dependerá del número de complejos PDB en el archivo de entrada.
- Cada fila del dataset contiene un identificador PDB y un vector de 100 elementos que representa las frecuencias n-hot.
"""


import pandas as pd

# Cargar los datos desde el archivo
file_path = './data/processed/fgs_pdb.csv'
output_path = './data/processed/fgs_n_hot.csv'
data = pd.read_csv(file_path)

# Contar la frecuencia de cada grupo funcional en la columna 'pseudo_smiles'
group_counts = data['pseudo_smiles'].value_counts()

# Seleccionar los 100 grupos funcionales más frecuentes
top_100_groups = group_counts.head(100).index

# Guardar los top 100 grupos funcionales en un archivo CSV
top_100_groups_df = pd.DataFrame({'pseudo_smiles': top_100_groups})
top_100_groups_path = './data/processed/top_100_fgs_ligandos_originales.csv'
top_100_groups_df.to_csv(top_100_groups_path, index=False)

# Crear una columna con listas de grupos funcionales
data['functional_groups'] = data['pseudo_smiles'].apply(lambda x: [x])

# Agrupar por 'Complejos PDB_ID' y consolidar los grupos funcionales en una lista única
grouped_data = data.groupby('Complejos PDB_ID')['functional_groups'].sum().reset_index()

# Contar la ocurrencia de cada grupo funcional en las listas
grouped_data['n_hot_fgs'] = grouped_data['functional_groups'].apply(
    lambda x: [x.count(group) for group in top_100_groups]
)

# Guardar los resultados en un archivo CSV
grouped_data[['Complejos PDB_ID', 'n_hot_fgs']].to_csv(output_path, index=False)