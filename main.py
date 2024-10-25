# main.py
from scripts import interacciones_proteina
from scripts import descarga_estructuras_pdb
from scripts import extract_ligands_uniprot
from scripts import fgs

# Ejecutar el script interacciones_proteina para generar el dataset principal
interacciones_proteina.main_interacciones_proteina()

# Ejecutar el script descarga_estructuras_pdb para descargar las estructuras CIF de PDB
descarga_estructuras_pdb.main_descarga_estructuras()

# Ejecutar el script extract_ligands_uniprot para extraer IDs de UniProt y ligandos
extract_ligands_uniprot.main_extract_ligands_uniprot()

# Ejecutar el script fgs para generar los pseudo-SMILES y grupos funcionales
fgs.process_pseudo_smiles()

print("Ejecución completa de todos los scripts.")
