"""
##########################
### Documentación del Dataset Final
##########################

Descripción General:
El dataset final combina información de múltiples fuentes relacionadas con compuestos químicos, proteínas objetivo y sus actividades biológicas.
Este dataset se utiliza principalmente para tareas de modelado químico-biológico, como predicción de actividad molecular o diseño de fármacos.

Origen de los Datos:
1. `dataset_with_pchembl`:
   - Datos de actividad biológica (pChEMBL) de ligandos.
   - Incluye identificadores únicos de complejos y proteínas (PDB_ID y UNIPROT_ID).
2. `one_hot_binding`:
   - Representación en formato one-hot encoding de características del sitio de unión de proteínas.
   - Procesado desde `binding_site_one_hot.csv`.
3. `n_hot_fgs`:
   - Representación funcional de grupos químicos (*functional groups*) en los complejos PDB.
   - Derivado de `fgs_pdb_n_hot_fgs.csv`.
4. `ligand_decoy`:
   - Información sobre ligandos señuelo (*decoys*) y sus características funcionales.
   - Derivado de `fgs_decoy_n_hot.csv`.

Objetivo del Dataset:
- Facilitar el desarrollo de modelos predictivos que relacionen características químicas y estructurales con la actividad biológica.
- Proporcionar un conjunto consolidado para estudios comparativos entre ligandos activos y decoys.
- Permitir el análisis de relaciones entre proteínas y ligandos mediante representación binaria y valores de actividad.

Estructura del Dataset Final:
- `pChEMBL`: Valor medio de pChEMBL para ligandos activos (0 para decoys).
- `Actividad`: Indicador binario de actividad (1 para ligandos activos, 0 para decoys).
- `n_hot_fgs`: Representación n-hot de los grupos funcionales asociados a los complejos.
- `One_Hot_Encoding`: Representación one-hot de las características del sitio de unión.

Tamaño del Dataset:
- Número total de filas: Depende de la combinación de ligandos activos y decoys seleccionados.
- Columnas principales: `pChEMBL`, `Actividad`, `n_hot_fgs`, `One_Hot_Encoding`.

Ubicación del Archivo Generado:
- El dataset final es guardado como un archivo CSV en la ruta especificada en `output_path`.
"""

# Importar librerías necesarias
import pandas as pd
from rdkit import Chem
from chembl_structure_pipeline import standardizer as sdz

# Función para convertir InChI a SMILES
def inchi_to_smiles(inchi):
    try:
        mol = Chem.MolFromInchi(inchi)
        return Chem.MolToSmiles(mol) if mol else None
    except Exception:
        return None

# Función para estandarizar SMILES
def standardize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            standardized_mol = sdz.standardize_mol(mol)
            return Chem.MolToSmiles(standardized_mol, isomericSmiles=False)
        return None
    except Exception:
        return None
    
# Cargar datasets
dataset_with_pchembl = pd.read_csv("./data/processed/dataset_with_pchembl.csv")
one_hot_binding = pd.read_csv("./data/processed/binding_site_one_hot.csv")
n_hot_fgs = pd.read_csv("./data/processed/fgs_pdb_n_hot_fgs.csv")
ligand_decoy = pd.read_csv("./data/processed/fgs_decoy_n_hot.csv")
output_path = "./data/processed/dataset_final.csv"

# Procesar dataset_with_pchembl
dataset_with_pchembl['Ligand_SMILES'] = dataset_with_pchembl['InChI'].apply(inchi_to_smiles)
dataset_with_pchembl['Actividad'] = 1
dataset_with_pchembl = dataset_with_pchembl[["Complejos PDB_ID", "UNIPROT_ID", "Ligand_SMILES", "pChEMBL", "Actividad"]]

#Procesar one_hot_binding
one_hot_binding.rename(columns={"PDB_ID": "Complejos PDB_ID"}, inplace=True)
one_hot_binding.drop(columns=["Binding_Site_Sequence_Ordered"], inplace=True)

# Fusionar datasets
merged_dataset = (
    dataset_with_pchembl
    .merge(n_hot_fgs, on='Complejos PDB_ID', how='inner')
    .merge(one_hot_binding, on='Complejos PDB_ID', how='inner')
)

# Estandarizar SMILES
for dataset in [merged_dataset, ligand_decoy]:
    dataset['Standardized_Ligand_SMILES'] = dataset['Ligand_SMILES'].apply(standardize_smiles)

# Filtrar SMILES correlacionados
correlated_smiles = set(ligand_decoy['Standardized_Ligand_SMILES']).intersection(
    merged_dataset['Standardized_Ligand_SMILES']
)
ligand_decoy = ligand_decoy[ligand_decoy['Standardized_Ligand_SMILES'].isin(correlated_smiles)]

# Agrupar datos de merged_dataset
grouped_dataset = merged_dataset.groupby(
    ['UNIPROT_ID', 'Standardized_Ligand_SMILES', 'One_Hot_Encoding']
).agg({
    'Complejos PDB_ID': lambda x: list(set(x)),
    'pChEMBL': 'mean'
}).reset_index()

# Mapear valores al dataset de decoys
ligand_decoy = ligand_decoy.merge(grouped_dataset, on='Standardized_Ligand_SMILES', how='left')

# Asignar valores por defecto
ligand_decoy['pChEMBL'] = 0
ligand_decoy['Actividad'] = 0

# Reducir el dataset a una muestra representativa
ligand_decoy_reduced = ligand_decoy.sample(frac=0.3, random_state=42)

# Seleccionar columnas finales para ligand_decoy_reduced
ligand_decoy_reduced = ligand_decoy_reduced[[
    'Complejos PDB_ID', 'UNIPROT_ID', 'Ligand_SMILES', 'pChEMBL',
    'Actividad', 'n_hot_fgs', 'One_Hot_Encoding', 'Standardized_Ligand_SMILES'
]]

# Concatenar datasets
dataset_final = pd.concat([merged_dataset, ligand_decoy_reduced], ignore_index=True)[
    ["pChEMBL", "Actividad", "n_hot_fgs", "One_Hot_Encoding"]
]

# Guardar el dataset final
dataset_final.to_csv(output_path, index=False)