import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

from utils import utils


if __name__ == "__main__":
    species = ['co', '12c16o']

    # import data
    df_fc = pd.read_csv("./data/abundance_colliders.csv", index_col="collider")

    df_collisions = utils.read_collisions(species)
    df_collisions = utils.process_collisions(df_collisions)

    df_lines = utils.read_lines(species)


    nu_list = list(df_lines['nu'].unique())
    kind_list = list(df_fc.columns)

    utils.plot_critical_density(species, nu_list[0], kind_list[1:], df_fc, df_collisions, df_lines)