import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import seaborn as sns

# Función para calcular estadísticas generales de fragmentos funcionales (FGs)
def calcular_estadisticas_fgs(fgs):
    # Utilizar la primera columna como referencia
    mol_col = fgs.iloc[:, 0]
    # Contar el número de moléculas únicas
    n_mols = len(mol_col.unique())
    # Crear un contador de FGs para contar ocurrencias
    fg_counter = Counter(fgs['pseudo_smiles'])
    # Contar el número total de fragmentos
    n_total_fgs = len(fgs)
    # Calcular el promedio de FGs por molécula
    fg_per_mol = n_total_fgs / n_mols

    # Contar FGs únicos en cada molécula (presencia/ausencia)
    binary_fg = fgs.groupby(mol_col)['pseudo_smiles'].nunique().sum()
    fg_bin_per_mol = binary_fg / n_mols

    # Calcular el número de FGs únicos en todo el dataset
    unique_fgs = len(fg_counter)
    unique_fg_per_mol = unique_fgs / n_mols

    # Crear una tabla con las estadísticas calculadas
    stats_table1 = {
        'n': [n_mols],
        'Total FGs': [n_total_fgs],
        'FGs por mol': [fg_per_mol],
        'FGs binarios': [binary_fg],
        'FGs binarios por mol': [fg_bin_per_mol],
        'FGs únicos': [unique_fgs],
        'FGs únicos por mol': [unique_fg_per_mol]
    }

    # Convertir a DataFrame para facilitar la visualización y uso posterior
    stats_table1_df = pd.DataFrame(stats_table1)

    return stats_table1_df

# Función para graficar los 20 fragmentos funcionales principales
def plot_top_20_fgs(fgs):
    # Crear un contador de FGs para contar ocurrencias
    fg_counter = Counter(fgs['pseudo_smiles'])
    # Obtener los 20 fragmentos funcionales más comunes
    top_20_fgs = fg_counter.most_common(20)
    fg_names, fg_counts = zip(*top_20_fgs)  # Desempaquetar nombres y conteos

    # Crear el gráfico de barras
    plt.figure(figsize=(10, 6))
    bars = plt.bar(fg_names, fg_counts)
    plt.xlabel('Fragmentos Funcionales (FG)')
    plt.ylabel('Número de apariciones')
    plt.title('Distribución de los 20 principales FGs')
    plt.xticks(rotation=90)

    # Añadir etiquetas con el conteo sobre cada barra
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), va='bottom')

    plt.tight_layout()
    plt.show()

# Función para calcular estadísticas de fragmentos aromáticos y con heteroátomos
def calcular_estadisticas_aromaticos_heteroatomos(fgs):
    # Definir una función para identificar si un FG es aromático
    def is_aromatic(smiles):
        return 'ar' in smiles

    # Definir una función para identificar si un FG contiene heteroátomos
    def is_heteroatom(smiles):
        heteroatoms = ['O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
        return any(atom in smiles for atom in heteroatoms)

    # Crear columnas en el DataFrame para marcar fragmentos aromáticos y con heteroátomos
    fgs['Aromatic_FG'] = fgs['pseudo_smiles'].apply(is_aromatic)
    fgs['Heteroatom_FG'] = fgs['pseudo_smiles'].apply(is_heteroatom)

    # Utilizar la primera columna como referencia
    mol_col = fgs.iloc[:, 0]
    # Agrupar por molécula (primera columna) y contar fragmentos aromáticos y con heteroátomos
    grouped = fgs.groupby(mol_col)
    aromatic_fgs_bin = grouped['Aromatic_FG'].sum()
    heteroatom_fgs_bin = grouped['Heteroatom_FG'].sum()

    # Calcular fragmentos únicos de aromáticos y con heteroátomos por molécula
    unique_aromatic_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'Aromatic_FG']].nunique())
    unique_heteroatom_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'Heteroatom_FG']].nunique())

    # Calcular valores promedio por molécula
    n_mols = len(mol_col.unique())
    ar_fgs_bin_per_mol = aromatic_fgs_bin.sum() / n_mols
    het_fgs_bin_per_mol = heteroatom_fgs_bin.sum() / n_mols
    unique_ar_fgs_per_mol = unique_aromatic_fgs.sum() / n_mols
    unique_het_fgs_per_mol = unique_heteroatom_fgs.sum() / n_mols

    # Crear tabla de estadísticas para los fragmentos aromáticos y con heteroátomos
    stats_table2 = {
        'Ar FGs (bin) / mol': [ar_fgs_bin_per_mol],
        '# Ar FGs (un)': [unique_aromatic_fgs.sum()],
        'Ar FGs (un) / mol': [unique_ar_fgs_per_mol],
        'Het FGs (bin) / mol': [het_fgs_bin_per_mol],
        '# Het FGs (un)': [unique_heteroatom_fgs.sum()],
        'Het FGs (un) / mol': [unique_het_fgs_per_mol]
    }

    stats_table2_df = pd.DataFrame(stats_table2)

    return stats_table2_df

def calcular_estadisticas_heteroatomos(fgs):
    # Funciones para identificar los FGs que contienen O, N, S, P y halógenos (X)
    def contains_oxygen(smiles):
        return 'O' in smiles

    def contains_nitrogen(smiles):
        return 'N' in smiles

    def contains_sulfur(smiles):
        return 'S' in smiles

    def contains_phosphorus(smiles):
        return 'P' in smiles

    def contains_halogen(smiles):
        # Halógenos comunes: F, Cl, Br, I
        return any(halogen in smiles for halogen in ['F', 'Cl', 'Br', 'I'])

    # Utilizar la primera columna como referencia
    mol_col = fgs.iloc[:, 0]

    # Crear las columnas en el dataframe para cada tipo de FG
    fgs['O_FG'] = fgs['pseudo_smiles'].apply(contains_oxygen)
    fgs['N_FG'] = fgs['pseudo_smiles'].apply(contains_nitrogen)
    fgs['S_FG'] = fgs['pseudo_smiles'].apply(contains_sulfur)
    fgs['P_FG'] = fgs['pseudo_smiles'].apply(contains_phosphorus)
    fgs['X_FG'] = fgs['pseudo_smiles'].apply(contains_halogen)

    # Agrupar por molécula (PDB_ID)
    grouped = fgs.groupby(mol_col)

    # Calcular los conteos binarios y únicos
    o_fgs_bin = grouped['O_FG'].sum()  # Conteo binario de FGs que contienen O
    n_fgs_bin = grouped['N_FG'].sum()  # Conteo binario de FGs que contienen N
    s_fgs_bin = grouped['S_FG'].sum()  # Conteo binario de FGs que contienen S
    p_fgs_bin = grouped['P_FG'].sum()  # Conteo binario de FGs que contienen P
    x_fgs_bin = grouped['X_FG'].sum()  # Conteo binario de FGs que contienen halógenos

    # Calcular los fragmentos únicos
    unique_o_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'O_FG']].nunique())
    unique_n_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'N_FG']].nunique())
    unique_s_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'S_FG']].nunique())
    unique_p_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'P_FG']].nunique())
    unique_x_fgs = grouped['pseudo_smiles'].apply(lambda x: x[fgs.loc[x.index, 'X_FG']].nunique())

    # Número total de moléculas
    n_mols = len(mol_col.unique())

    # Calcular los valores por molécula y las fracciones
    o_fgs_bin_per_mol = o_fgs_bin.sum() / n_mols
    n_fgs_bin_per_mol = n_fgs_bin.sum() / n_mols
    s_fgs_bin_per_mol = s_fgs_bin.sum() / n_mols
    p_fgs_bin_per_mol = p_fgs_bin.sum() / n_mols
    x_fgs_bin_per_mol = x_fgs_bin.sum() / n_mols

    frac_o_fgs_un = unique_o_fgs.sum() / n_mols
    frac_n_fgs_un = unique_n_fgs.sum() / n_mols
    frac_s_fgs_un = unique_s_fgs.sum() / n_mols
    frac_p_fgs_un = unique_p_fgs.sum() / n_mols
    frac_x_fgs_un = unique_x_fgs.sum() / n_mols

    # Crear la tabla de estadísticas
    stats_table3 = pd.DataFrame({
        'O FGs (bin) / mol': [o_fgs_bin_per_mol],
        'Frac O FGs (un)': [frac_o_fgs_un],
        'N FGs (bin) / mol': [n_fgs_bin_per_mol],
        'Frac N FGs (un)': [frac_n_fgs_un],
        'S FGs (bin) / mol': [s_fgs_bin_per_mol],
        'Frac S FGs (un)': [frac_s_fgs_un],
        'P FGs (bin) / mol': [p_fgs_bin_per_mol],
        'Frac P FGs (un)': [frac_p_fgs_un],
        'X FGs (bin) / mol': [x_fgs_bin_per_mol],
        'Frac X FGs (un)': [frac_x_fgs_un]
    })

    return stats_table3

# Función para visualizar la co-ocurrencia de los 20 fragmentos funcionales más comunes
def visualizar_coocurrencia_top_20_fgs(fgs):
    # Utilizar la primera columna como referencia
    mol_col = fgs.iloc[:, 0]

    # Crear una tabla de presencia/ausencia de FGs por molécula usando la primera columna
    fg_presence = pd.crosstab(mol_col, fgs['pseudo_smiles'])

    # Calcular la matriz de co-ocurrencia usando producto matricial
    co_occurrence_matrix = np.dot(fg_presence.T, fg_presence)

    # Convertir la matriz a DataFrame para etiquetar correctamente los fragmentos funcionales
    co_occurrence_df = pd.DataFrame(co_occurrence_matrix, index=fg_presence.columns, columns=fg_presence.columns)

    # Seleccionar los 20 fragmentos más comunes y filtrar la matriz de co-ocurrencia
    top_20_fgs = fgs['pseudo_smiles'].value_counts().index[:20]
    co_occurrence_top_20 = co_occurrence_df.loc[top_20_fgs, top_20_fgs]

    # Convertir la matriz a float para permitir NaN
    co_occurrence_top_20 = co_occurrence_top_20.astype(float)

    # Eliminar la diagonal (asignarla a np.nan para que no se considere en la escala de color)
    np.fill_diagonal(co_occurrence_top_20.values, np.nan)

    # Visualizar la matriz de co-ocurrencia con un mapa de calor
    plt.figure(figsize=(10, 8))
    sns.heatmap(co_occurrence_top_20, cmap='YlGnBu', annot=False)
    plt.title("Matriz de Co-ocurrencia de los Top 20 Grupos Funcionales")
    plt.xlabel("Grupos Funcionales")
    plt.ylabel("Grupos Funcionales")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
