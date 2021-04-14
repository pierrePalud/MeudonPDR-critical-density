import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

import glob
import re
import os

def read_lines(species):
    """reads the file that contains the Einstein coefficients
    of the radiative de-excitation

    Parameters
    ----------
    species : list or str
        name or list of names of the one species to consider
        (if list of names : only one of them can correspond to a file
        in ./data/Lines) 

    Returns
    -------
    df_lines : pandas.DataFrame
        content of the table in the file indicated by the species name
    """
    if isinstance(species, str):
        species = [species]
    assert isinstance(species, list)

    list_filenames = []
    for species_name in species:
        list_filenames += sorted(glob.glob(f"./data/Lines/line_{species_name}.dat"))

    if len(list_filenames) == 0:
        raise ValueError(f"There are no files in ./data/Lines/ that contain one of the elements {species}")

    if len(list_filenames) > 1:
        msg = f"There are more than 1 file in ./data/Lines/ that contain one of the elements {species} \n"
        msg += "Please only indicate one chemical species."
        raise ValueError(msg)


    list_lines = []
    for filename in list_filenames:
        # reads in the file the number of rows in the table (first line)
        with open(filename) as f:
            first_line =  f.readline()
            num_rows = int(first_line.split('#')[1].replace(' ', ''))

        df_lines = pd.read_csv(
            filename, skiprows=2, sep=r"\s+", engine="python", comment='#--', nrows=num_rows
        )
        df_lines.columns = list(df_lines.columns[1:]) + ["unit"]

        # withdraw the semicolumn in quantum numbers columns
        for col in df_lines.select_dtypes(include="object").columns:
            df_lines[col] = df_lines[col].apply(lambda x: x.replace(";", ""))

            try:
                df_lines[col] = df_lines[col].astype(int)
            except:
                pass

        # indicate the file from which the data is from
        df_lines['file'] = filename

        list_lines.append(df_lines)

    # group results in one dataframe
    df_lines = pd.concat(list_lines)    
    return df_lines


def compute_sum_einstein_coefs(df_lines, nu):
    """sums the Einstein coefficients of all transitions from 
    a fixed upper level to lower levels 

    Parameters
    ----------
    df_lines : pandas.DataFrame
        the content of the ./data/Lines file corresponding to the desired species

    nu : int
        index of the upper level

    Returns
    -------
    sum_einstein_coefs : float
        sum of the Einstein coefficients of all the transitions from the nu 
        level to lower levels 
    """
    # get the name of the column with einstein coefficient of transition
    # (this column is sometimes called "Aij(s-1)" or "Aein(s-1)")
    coef_einstein_colname = re.search(r"A[a-z]+\(s-1\)", ''.join(list(df_lines.columns)))[0]

    # process lines data (sum over nl<nu)
    df_sum_einstein_coefs = df_lines[["nu", coef_einstein_colname]]
    df_sum_einstein_coefs = df_sum_einstein_coefs[df_sum_einstein_coefs["nu"] == nu]
    df_sum_einstein_coefs = df_sum_einstein_coefs.groupby('nu')[coef_einstein_colname].sum()

    # extract result 
    sum_einstein_coefs = df_sum_einstein_coefs.at[nu]
    return sum_einstein_coefs



def read_collisions(species):
    if isinstance(species, str):
        species = [species]
    assert isinstance(species, list)

    # get all collisions files in which the spcecies appear
    list_filenames = []
    for species_name in species:
        list_filenames += sorted(glob.glob(f"./data/Collisions/*{species_name}[_.]dat"))

    list_collisions = []
    for filename in list_filenames:
        # define header
        idx_header = 8
        with open(filename) as f:
            for i in range(idx_header - 1):
                val = f.readline()

                # number of rows to read / size of the table in 4th line
                if i==3:
                    nrows = int(val)

            header = f.readline()
            header = re.split(r'\s+', header) # split
            header = [float(col) for col in header if len(col) > 0]
            header = ["nu", "nl"] + header


        # read content
        df_collisions = pd.read_csv(
            filename, skiprows=9, sep=r"\s+", engine="python", comment='#-',
            header=None, nrows=nrows,
            skip_blank_lines=True
        )
        if len(header) + 1 == len(df_collisions.columns):
            df_collisions = df_collisions.drop([0], 1)

        df_collisions.columns = header

        df_collisions['collider'] = filename.split("_")[1]
        df_collisions = df_collisions.set_index(["nu", "nl", "collider"])

        list_collisions.append(df_collisions)

    df_collisions = pd.concat(list_collisions)
    df_collisions = df_collisions.reindex(sorted(df_collisions.columns), axis=1)
    return df_collisions


def process_collisions(df_collisions):

    # linear interpolation collision coefficients
    # (only works if there are known values for lower AND higher temperatures for that collider)
    # (otherwise the whole column is dropped)
    df_collisions = df_collisions.interpolate(axis=1)
    df_collisions = df_collisions.dropna(axis=1)

    return df_collisions



def compute_critical_density(kind, nu, sum_einstein_coefs, df_collisions, df_fc):
    # process collisions (for each collisioner, sum over nl<nu)
    # according to the chosen kind of medium
    df_collisions_per_collisioner = df_collisions.loc[nu]
    df_collisions_per_collisioner = df_collisions_per_collisioner.groupby(level=[1]).sum() # groupby collider


    df_collisions_global = pd.merge(
        df_collisions_per_collisioner, 
        df_fc[[kind]], 
        left_index=True, right_index=True
    )

    # conclude computation of denominator
    for col in df_collisions_global.columns:
        if col != kind:
            df_collisions_global[col] *= df_collisions_global[kind] 

    df_collisions_global = df_collisions_global.drop(kind, 1)
    df_collisions_global = df_collisions_global.sum() 

    # create divide the numerator by the denominator
    df_ncr = sum_einstein_coefs / df_collisions_global 
    df_ncr = df_ncr.apply(lambda x: np.nan if np.isinf(x) else x)
    # for col in df_ncr.columns: 
    #     df_ncr[col] = sum_einstein_coefs / df_collisions_global[col]
        # df_ncr[col] = df_ncr[col].apply(lambda x: np.nan if np.isinf(x) else x)


    return df_ncr



def plot_critical_density(species, nu, kinds, df_fc, df_collisions, df_lines):
    fig, ax = plt.subplots(figsize=(8,6))

    # compute the critical density for a set of temperatures and plot them
    sum_einstein_coefs = compute_sum_einstein_coefs(df_lines, nu)
    for kind in kinds:
        df_ncr = compute_critical_density(kind, nu, sum_einstein_coefs, df_collisions, df_fc)
        plt.loglog(df_ncr.index, df_ncr.values, "+-", label=kind)

    # set title : get quantum numbers of level u
    quantum_numbers_colnames = list(df_lines.columns)
    i_start = quantum_numbers_colnames.index('quant:')
    i_end = quantum_numbers_colnames.index('info:')
    quantum_numbers_colnames = quantum_numbers_colnames[i_start+1:i_end]
    quantum_numbers_colnames = [col for col in quantum_numbers_colnames if "u" in col]

    df_quantum_numbers_of_level = df_lines.set_index("nu")[quantum_numbers_colnames].drop_duplicates()
    dict_quantum_numbers_of_level = df_quantum_numbers_of_level.loc[nu].to_dict()

    plt.title(f"critical densities of {species[0]} for nu = {nu} ({dict_quantum_numbers_of_level})")

    # conclude figure layout
    plt.xlabel('Temperature (K)')
    plt.ylabel("density (cm-3)")
    plt.legend()
    plt.grid()
    plt.show()