import os
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

# Definir el mapa de aminoácidos
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
amino_to_index = {aa: idx for idx, aa in enumerate(amino_acids)}

# Función para convertir una secuencia en one-hot encoding
def sequence_to_one_hot(sequence, amino_to_index):
    one_hot = np.zeros((len(sequence), len(amino_to_index)), dtype=int)
    for i, aa in enumerate(sequence):
        if aa in amino_to_index:
            one_hot[i, amino_to_index[aa]] = 1
    return one_hot

# Función para obtener el sitio de unión de la primera cadena con el ligando
def get_binding_site_first_chain(cif_file, ligand_name, radius=4):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", cif_file)
    except Exception as e:
        print(f"Error procesando {cif_file}: {e}")
        return None

    model = structure[0]

    # Buscar el ligando en la primera cadena que lo contenga
    ligand_atoms = None
    chain_id_with_ligand = None

    for chain in model:
        for residue in chain:
            if residue.resname == ligand_name:
                ligand_atoms = list(residue.get_atoms())
                chain_id_with_ligand = chain.id
                break
        if ligand_atoms:
            break

    if not ligand_atoms:
        print(f"Ligando {ligand_name} no encontrado en {cif_file}")
        return None

    # Buscar residuos cercanos al ligando en la cadena seleccionada
    ns = NeighborSearch(list(model.get_atoms()))
    close_atoms = []
    for atom in ligand_atoms:
        close_atoms.extend(ns.search(atom.coord, radius))

    # Contar átomos por residuo
    residue_atom_counts = {}
    for close_atom in close_atoms:
        parent_residue = close_atom.get_parent()
        if is_aa(parent_residue, standard=True):
            if parent_residue not in residue_atom_counts:
                residue_atom_counts[parent_residue] = 0
            residue_atom_counts[parent_residue] += 1

    # Filtrar residuos con al menos 2 átomos cercanos
    filtered_residues = [residue for residue, count in residue_atom_counts.items() if count >= 2]

    # Generar lista de residuos en formato "resname+id" y ordenarlos
    binding_site = list(set(f"{res.resname}{res.id[1]}" for res in filtered_residues))
    binding_site_ordered = sorted(binding_site, key=lambda x: int(''.join(filter(str.isdigit, x))))

    # Convertir a secuencia de aminoácidos
    amino_acid_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    sequence_ordered = ''.join([amino_acid_map.get(res[:3], 'X') for res in binding_site_ordered])

    return {
        "Chain_ID": chain_id_with_ligand,
        "Binding_Site_Sequence_Ordered": sequence_ordered
    }

# Función para dividir un DataFrame en lotes
def split_dataframe(df, batch_size):
    """
    Divide un DataFrame en lotes de tamaño batch_size.
    """
    for i in range(0, len(df), batch_size):
        yield df.iloc[i:i + batch_size]

# Ruta de entrada y salida
file_path = "./data/processed/filtered_extract_ligands_uniprot.csv"
cif_folder = "./data/estructuras_pdb"
output_file_one_hot = "./data/processed/binding_site_one_hot.csv"

# Tamaño de cada lote
batch_size = 100  # Ajusta según tus necesidades

# Cargar el archivo CSV
proteins = pd.read_csv(file_path)

# Eliminar duplicados basados en la columna 'InChI' y 'UNIPROT_ID'
proteins = proteins.drop_duplicates(subset=['InChI', 'UNIPROT_ID'])

# Eliminar filas con 'Ligand PDB_ID' nulo
proteins = proteins.dropna(subset=['Ligand PDB_ID'])

# Crear carpeta de salida si no existe
os.makedirs(os.path.dirname(output_file_one_hot), exist_ok=True)

# Acumular resultados de todos los lotes
all_results = []

# Procesar proteínas por lotes
batch_counter = 1
for batch in split_dataframe(proteins, batch_size):
    print(f"Procesando lote {batch_counter}...")

    for _, row in batch.iterrows():
        pdb_id = row['Complejos PDB_ID']
        ligand_id = row['Ligand PDB_ID']

        if pd.isnull(ligand_id) or not isinstance(ligand_id, str):
            continue

        cif_file = os.path.join(cif_folder, f"{pdb_id}.cif")
        if not os.path.isfile(cif_file):
            continue

        ligand_name = ligand_id[:3]
        result = get_binding_site_first_chain(cif_file, ligand_name)

        if result and result["Binding_Site_Sequence_Ordered"]:
            one_hot_vector = sequence_to_one_hot(result["Binding_Site_Sequence_Ordered"], amino_to_index)
            all_results.append({
                "PDB_ID": pdb_id,
                "Binding_Site_Sequence_Ordered": result["Binding_Site_Sequence_Ordered"],
                "One_Hot_Encoding": one_hot_vector.tolist()
            })

    print(f"Lote {batch_counter} procesado.")
    batch_counter += 1

# Guardar todos los resultados en un único archivo CSV
all_results_df = pd.DataFrame(all_results)
all_results_df.to_csv(output_file_one_hot, index=False)

print(f"Resultados de one-hot encoding guardados en: {output_file_one_hot}")
