# -*- coding: utf-8 -*-
"""
Script: descarga_estructuras_pdb.py
@author: Gabriel

Descripción:
Este script descarga estructuras de proteínas en formato CIF desde el RCSB Protein Data Bank (PDB),
utilizando un archivo de texto que contiene una lista de identificadores PDB.
El flujo incluye los siguientes pasos:

1. Leer los identificadores PDB desde un archivo de texto.
2. Descargar cada archivo CIF correspondiente a los PDB IDs desde el servidor de RCSB PDB.
3. Guardar los archivos descargados en una carpeta específica.
4. Registrar en un archivo de log los archivos descargados y aquellos que no se pudieron descargar.

Salidas:
- Archivos CIF descargados en la carpeta `data/estructuras_pdb/`.
- Archivos de log con detalles sobre los archivos descargados correctamente y aquellos que presentaron errores.
"""

import os
import requests

# Función para leer los identificadores PDB desde un archivo de texto
def leer_pdb_ids(file_path):
    """
    Lee un archivo de texto que contiene una lista de PDB IDs y los retorna en una lista.

    Parámetros:
        file_path (str): Ruta al archivo de texto con los PDB IDs.

    Retorna:
        list: Lista de PDB IDs leídos del archivo, sin espacios en blanco.
    """
    with open(file_path, 'r') as file:
        pdb_ids = [pdb_id.strip() for pdb_id in file.read().split('\n') if pdb_id.strip()]
    return pdb_ids

# Función para descargar archivos CIF desde el RCSB PDB
def descargar_estructuras_pdb(pdb_ids, path2save, timeout=30):
    """
    Descarga archivos CIF de proteínas utilizando los PDB IDs, desde la API de RCSB PDB.

    Parámetros:
        pdb_ids (list): Lista de PDB IDs que se desean descargar.
        path2save (str): Ruta donde se guardarán los archivos descargados.
        timeout (int): Tiempo de espera para la solicitud en segundos (default: 30 segundos).

    Retorna:
        tuple: Número de archivos descargados correctamente y una lista de PDB IDs que fallaron.
    """
    api_url = 'https://files.rcsb.org/download'
    
    # Crear el directorio de descarga si no existe
    if not os.path.exists(path2save):
        os.makedirs(path2save)

    downloaded_count = 0  # Contador de archivos descargados exitosamente
    errores = []  # Lista de PDB IDs que no se pudieron descargar

    for pdb_id in pdb_ids:
        # Definir URL de descarga y ruta de archivo de destino
        url = f"{api_url}/{pdb_id}.cif"
        file_path = os.path.join(path2save, f"{pdb_id}.cif")

        # Solo descargar si el archivo no existe en la carpeta
        if not os.path.exists(file_path):
            try:
                # Realizar la solicitud de descarga
                response = requests.get(url, timeout=timeout)
                response.raise_for_status()  # Lanza error si la respuesta es un error HTTP
                with open(file_path, "wb") as cif_file:
                    cif_file.write(response.content)  # Guardar el archivo CIF descargado
                print(f"Downloaded {pdb_id}.cif")
                downloaded_count += 1
            except requests.exceptions.RequestException as e:
                # Manejar y registrar errores de descarga
                print(f"Failed to download {pdb_id}. Error: {e}")
                errores.append(pdb_id)
        else:
            print(f"{pdb_id}.cif already exists.")  # Notificación si el archivo ya existe
    
    print(f"Total number of files downloaded: {downloaded_count}")
    return downloaded_count, errores

# Función para guardar los errores de descarga en un archivo
def guardar_errores(errores, error_file_path):
    """
    Guarda en un archivo de texto los PDB IDs que no se pudieron descargar.

    Parámetros:
        errores (list): Lista de PDB IDs que no se pudieron descargar.
        error_file_path (str): Ruta del archivo donde se guardarán los PDB IDs con errores.
    """
    with open(error_file_path, 'w') as error_file:
        if errores:
            error_file.write("Los siguientes complejos no se pudieron descargar correctamente:\n")
            for nombre in errores:
                error_file.write(f"{nombre}\n")
        else:
            error_file.write("Todos los complejos se descargaron correctamente.\n")

# Función principal que coordina el proceso completo de descarga
def main_descarga_estructuras():
    """
    Ejecuta el flujo completo de descarga de estructuras PDB:
    1. Lee los PDB IDs desde un archivo de texto.
    2. Descarga los archivos CIF de cada ID.
    3. Registra en archivos de log los resultados de las descargas.
    """
    # Definir rutas de archivos de entrada y salida
    pdb_ids_file_path = "./data/processed/complejos_PDB.txt"  # Archivo con los IDs PDB a descargar
    path2save = "./data/estructuras_pdb/"                     # Carpeta de destino para los archivos CIF descargados
    log_file_path = "./data/processed/salida_terminal.txt"    # Archivo de log para resumen de la descarga
    error_file_path = "./data/processed/errores.txt"          # Archivo de log para los PDB IDs con errores
    
    # Leer los PDB IDs desde el archivo
    print("Leyendo PDB IDs...")
    pdb_ids = leer_pdb_ids(pdb_ids_file_path)
    
    # Iniciar la descarga de estructuras PDB
    print("Iniciando la descarga de estructuras PDB...")
    descargados, errores = descargar_estructuras_pdb(pdb_ids, path2save, timeout=180)
    
    # Guardar los resultados de la descarga y los errores en archivos de log
    print("Guardando resultados en archivo de log...")
    guardar_errores(errores, error_file_path)
    
    # Guardar un resumen en el archivo de log
    with open(log_file_path, 'w') as log_file:
        log_file.write(f"Número total de proteínas a descargar: {len(pdb_ids)}\n")
        log_file.write(f"Número de archivos descargados correctamente: {descargados}\n")
        log_file.write(f"Número de archivos que no se pudieron descargar: {len(errores)}\n")

    print(f"Proceso completado. Log guardado en '{log_file_path}' y errores en '{error_file_path}'.")
