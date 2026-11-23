#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:52:25 2026

@author: mckayapugh
"""
from epydemix import EpiModel
from epydemix.visualization import plot_quantiles

from epydemix import load_predefined_model
from epydemix.utils import convert_to_2Darray, compute_simulation_dates
from epydemix.population import load_epydemix_population

# population = load_epydemix_population(population_name="South_Korea") # Could eventually use for age stratification, haven't tried it out yet

#### Paper's parameters 
"""
R = 1.46
Incubation period = 1.5 days
Infectious period = 3 days
Asx prop: 50% [35%, 65%]
Asx reduction ratio: 50% [35%, 65%]

VE = 60%
Immunity waning rate = 0.56% per day
N = 4,000 (2,000 per group)
"""

# !!! Missing non-trial compartments
# -- Model definition --
model = EpiModel(
    name="SEIR Model",
    compartments=["Sp","Ep", "Ip", "Ap","PI","R", # Placebo arm
                  "Svf", "Evf", "Ivf", "Avf", "VI",  # Vaccine arm: failure 
                  "V", "Sv", "Ev", "Iv", "Av"  # Vaccine arm: success
                  ],
    parameters={"beta": 0.486,
                "beta_a": 0.486*0.5,
                "gamma": 0.333,
                "sigma": 0.666,
                "w": 0.00
                }
)

# -- transition rules ---

# Placebo arm
# Transmissions: Infectious compartments
model.add_transition(source="Sp", target="Ep", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Sp", target="Ep", params=("beta_a", "Ap"), kind="mediated")
model.add_transition(source="Sp", target="Ep", params=("beta_a", "Avf"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta_a", "Av"), kind="mediated") 
# Latent transitions
model.add_transition(source="Ep", target="Ip", params="sigma", kind="spontaneous")
model.add_transition(source="Ep", target="Ap", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Ip", target="PI", params="gamma", kind="spontaneous") 
model.add_transition(source="Ap", target="R", params="gamma", kind="spontaneous")

# Vaccine arm: failure
# Transmissions: Infectious compartments
model.add_transition(source="Svf", target="Evf", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Svf", target="Evf", params=("beta_a", "Ap"), kind="mediated")
model.add_transition(source="Svf", target="Evf", params=("beta_a", "Avf"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta_a", "Av"), kind="mediated") 
# Latent transitions
model.add_transition(source="Evf", target="Ivf", params="sigma", kind="spontaneous")
model.add_transition(source="Evf", target="Avf", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Ivf", target="VI", params="gamma", kind="spontaneous") 
model.add_transition(source="Avf", target="R", params="gamma", kind="spontaneous")

# Vaccine arm: success
# Waning immunity
model.add_transition(source="V", target="Sv", params=("w"), kind="spontaneous") 
# Transmissions: Infectious compartments
model.add_transition(source="Sv", target="Ev", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Sv", target="Ev", params=("beta_a", "Ap"), kind="mediated")
model.add_transition(source="Sv", target="Ev", params=("beta_a", "Avf"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta_a", "Av"), kind="mediated") 
# Latent transitions
model.add_transition(source="Ev", target="Iv", params="sigma", kind="spontaneous")
model.add_transition(source="Ev", target="Av", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Iv", target="VI", params="gamma", kind="spontaneous") 
model.add_transition(source="Av", target="R", params="gamma", kind="spontaneous")



N = 100000
VE = 0.6
s0 = 0.99
i0 = 1-s0

# Run simulations
results = model.run_simulations(
    initial_conditions_dict={"Sp": s0*N, #"V": s0*1/2*VE*N, "Svf": s0*1/2*(1-VE)*N,
                             "Ap": i0*N},
    start_date="2024-01-01",
    end_date="2024-12-10",
    Nsim=100,
)

# Extract and plot quantiles of compartment counts
df_quantiles = results.get_quantiles_compartments()
plot_quantiles(df_quantiles, columns=[ "Sp_total", "Ep_total", "Ip_total","R_total"])
plot_quantiles(df_quantiles, columns=[ "Ip_total", "Ivf_total", "Iv_total"])

# -- Calculating VE estimate --
dfq = df_quantiles
finals = dfq[(dfq['date'] == "2024-09-10 00:00:00") & (dfq['quantile'] == 0.5)]

VEhat = 1 - finals['VI_total']/finals['PI_total']
print(VEhat)
