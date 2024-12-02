import pandas as pd

# Cargar los datos desde el archivo
file_path = './data/processed/fgs_pdb.csv'
output_path = './data/processed/fgs_n_hot.csv'
data = pd.read_csv(file_path)

# Contar la frecuencia de cada grupo funcional en la columna 'pseudo_smiles'
group_counts = data['pseudo_smiles'].value_counts()

# Seleccionar los 100 grupos funcionales más frecuentes
top_100_groups = group_counts.head(100).index

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