#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:11:39 2026

The way I programmed this is probably horrifying to someone in comp sci. Please
dont' judge too hard, I am but a lowly applied mathematician who messes around
with the code haphazardly until it works. Computational efficiency is none of 
my business.

Install from conda: conda install -c conda-forge epydemix

@author: mckayapugh
"""

from epydemix import EpiModel
from epydemix.visualization import plot_quantiles

from epydemix import load_predefined_model
from epydemix.utils import convert_to_2Darray, compute_simulation_dates
from epydemix.population import load_epydemix_population

import numpy as np


#### Paper's parameters 
"""
R = 1.46
Incubation period = 1.5 days
Infectious period = 3 days
Asx prop: 50% [35%, 65%]
Asx reduction ratio: 50% [35%, 65%]

VE = 60%
Immunity waning rate = 0.56% per day   
   --> average immunity duration 6mo 
N = 4,000 (2,000 per group)

Trial:
4k participants recruited over 4 weeks (3-4 per day)
"""

model = EpiModel(
    name="SEIR Model",
    compartments=["S", "E", "I", "A", "R", # general population
                  "Sp","Ep", "Ip", "PI", # Placebo arm
                  "Svf", "Evf", "Ivf", "VI",  # Vaccine arm: failure 
                  "V", "Sv", "Ev", "Iv"  # Vaccine arm: success
                  ],
    parameters={"beta": 0.486,
                "beta_a": 0.486*0.5,
                "gamma": 0.333,
                "sigma": 0.666,
                "w": 0.0056,
                "VE": 0.6
                }
)

# Including epydemix's age strat and contact matrix for South Korea. As is it has
# less compartments than the og paper and the contact matrix maybe looks a little 
# different, but I figure we can decide if that really matters later.
# Right now it doesn't really do anyhting bc I haven't set up age dependent transmission
population = load_epydemix_population(population_name="South_Korea")
model.set_population(population)

start = "2024-01-01"
end = "2024-03-01"
timeperiod = compute_simulation_dates(start_date=start, end_date=end)

rec_period = 7*4  # 4 weeks
recruitment = [1/3.5]*(rec_period) + [0]*(len(timeperiod)-rec_period)

model.add_parameter(parameter_name="recruitment_rate", value = recruitment)

### Add transitions: 

# -- General population
# Trial recruitment
model.add_transition(source="S", target="Sp", params="1/2*recruitment_rate", kind="spontaneous")
model.add_transition(source="S", target="Svf", params="1/2*(1-VE)*recruitment_rate", kind="spontaneous")
model.add_transition(source="S", target="V", params="1/2*VE*recruitment_rate", kind="spontaneous")

# Transmissions: Infectious compartments
model.add_transition(source="S", target="E", params=("beta", "I"), kind="mediated") 
model.add_transition(source="S", target="E", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="S", target="E", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="S", target="E", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="S", target="E", params=("beta_a", "A"), kind="mediated")
# Latent transitions
model.add_transition(source="E", target="I", params="sigma", kind="spontaneous")
model.add_transition(source="E", target="A", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="I", target="R", params="gamma", kind="spontaneous") 
model.add_transition(source="A", target="R", params="gamma", kind="spontaneous")

  
# -- Placebo arm
# Transmissions: Infectious compartments
model.add_transition(source="Sp", target="Ep", params=("beta", "I"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Sp", target="Ep", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Sp", target="Ep", params=("beta_a", "A"), kind="mediated")
# Latent transitions
model.add_transition(source="Ep", target="Ip", params="sigma", kind="spontaneous")
model.add_transition(source="Ep", target="A", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Ip", target="PI", params="gamma", kind="spontaneous") 

# -- Vaccine arm, failure
# Transmissions: Infectious compartments
model.add_transition(source="Svf", target="Evf", params=("beta", "I"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Svf", target="Evf", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Svf", target="Evf", params=("beta_a", "A"), kind="mediated")
# Latent transitions
model.add_transition(source="Evf", target="Ivf", params="sigma", kind="spontaneous")
model.add_transition(source="Evf", target="A", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Ivf", target="VI", params="gamma", kind="spontaneous") 

# -- Vaccine arm, success
# Waning immunity
model.add_transition(source="V", target="Sv", params=("w"), kind="spontaneous") 
# Transmissions: Infectious compartments
model.add_transition(source="Sv", target="Ev", params=("beta", "I"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta", "Ip"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta", "Ivf"), kind="mediated") 
model.add_transition(source="Sv", target="Ev", params=("beta", "Iv"), kind="mediated") 
# Transmissions: Asymptomatic compartments
model.add_transition(source="Sv", target="Ev", params=("beta_a", "A"), kind="mediated")
# Latent transitions
model.add_transition(source="Ev", target="Iv", params="sigma", kind="spontaneous")
model.add_transition(source="Ev", target="A", params="sigma", kind="spontaneous")
# Recoveries
model.add_transition(source="Iv", target="VI", params="gamma", kind="spontaneous")



N = 51738086  ## make sure this matches

s0 = 0.9999
i0 = 0.0001

# Run simulations
results = model.run_simulations(
    initial_conditions_dict={"S": s0*N, "I": i0*N},
    start_date= start,
    end_date= end,
    Nsim=100,
)

# Extract and plot quantiles of compartment counts
df_quantiles = results.get_quantiles_compartments()
plot_quantiles(df_quantiles, columns=[ "S_total", "E_total", "I_total","R_total"]) # Plot just gen pop
plot_quantiles(df_quantiles, columns=[ "Ip_total", "Ivf_total", "Iv_total"]) # Plot tial infection compartments

dfq = df_quantiles
finals = dfq[(dfq['date'] == end + "00:00:00") & (dfq['quantile'] == 0.5)]

VEhat = 1 - finals['VI_total']/finals['PI_total']
print(VEhat)
