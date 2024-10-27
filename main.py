# main.py
from scripts import interacciones_proteina
from scripts import descarga_estructuras_pdb
from scripts import extract_ligands_uniprot
from scripts.fgs import generate_fgs_pdb, generate_fgs_drugbank

# Ejecutar el script interacciones_proteina para generar el dataset principal
interacciones_proteina.main_interacciones_proteina()

# Ejecutar el script descarga_estructuras_pdb para descargar las estructuras CIF de PDB
descarga_estructuras_pdb.main_descarga_estructuras()

# Ejecutar el script extract_ligands_uniprot para extraer IDs de UniProt y ligandos
extract_ligands_uniprot.main_extract_ligands_uniprot()

# Ejecutar el script fgs para generar los pseudo-SMILES y grupos funcionales
# Generar los archivos fgs_pdb.csv y fgs_drugbank.csv usando la función general
generate_fgs_pdb("./data/processed/filtered_extract_ligands_uniprot.csv", "./data/processed/fgs_pdb.csv")
generate_fgs_drugbank("./data/processed/result_drugbank.csv", "./data/processed/fgs_drugbank.csv")

print("Ejecución completa de todos los scripts.")
