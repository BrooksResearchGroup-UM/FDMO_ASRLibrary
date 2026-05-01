import plotly.express as px
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## change current directory
os.chdir("/home/cdchiang/for_azamh/vae")

# Model directory
model_dir = "20221117_PF01494"

# Load model and the latent space representations of the training and testing data
with open("./{}/NoCV_latent_space_d2_layer4_w0.0005_b256_l0.0005_beta1_1000epoch_seed19.pkl".format(model_dir), 'rb') as file_handle:
    latent_space = pickle.load(file_handle)

# training datatset
key = latent_space['key'] 
mu = latent_space['mu']
sigma = latent_space['sigma']

# testing dataset
t_key = latent_space['t_key'] 
t_mu = latent_space['t_mu']
t_sigma = latent_space['t_sigma']

# Plot mapped wild type sequences in the latent space
fig = px.scatter(x=t_mu[:781,0], y=t_mu[:781,1], range_x = [-0.5, 1.5], 
                 range_y = [-3, -1], width=800, height=800)
fig.update_layout(
    title="PF01494 testing data (wild type)",
    xaxis_title="Z1",
    yaxis_title="Z2",
    #legend_title="Legend Title",
    font=dict(
        family="Arial, monospace",
        size=18,
        color="Black"))

fig.write_image("./{}/PF01494_testing_WT.png".format(model_dir))

# Plot mapped evolved sequences in the latent space
fig = px.scatter(x=t_mu[781:,0], y=t_mu[781:,1], range_x = [-0.5, 1.5], 
                 range_y = [-3, -1], width=800, height=800)
fig.update_layout(
    title="PF01494 testing data (evolved)",
    xaxis_title="Z1",
    yaxis_title="Z2",
    #legend_title="Legend Title",
    font=dict(
        family="Arial, monospace",
        size=18,
        color="Black"))

fig.write_image("./{}/PF01494_testing_evolved.png".format(model_dir))
