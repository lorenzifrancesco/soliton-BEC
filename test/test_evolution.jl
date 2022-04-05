using SolitonBEC
using Printf

### SIMULATE GENERATION OF SOLITON IN Li-7 CONDENSATE

## Simulation
std_sim = Simulation(300e-3, #T
  5e-4, #dt
  100e-6, #S
  1e-8,
  "GPE") #ds

npse_sim = Simulation(300e-3, #T
  5e-4, #dt
  100e-6, #S
  1e-8,
  "NPSE") #ds
  
## Apparata
std_apparatus = Apparatus(1e-26, #m
               1.4e-9, #as
               400, # omega_perp
               4e8) # N

## InitialStates      
InitialState1 =  InitialState(1,
                2e-12,
                0)    

## Configuration
configs = [(std_sim, std_apparatus, InitialState1), 
(npse_sim, std_apparatus, InitialState1)]
mem_limit = 15000000000 #byte

for (sim, app, state) in configs
  estimate = mem_estimate(sim)
  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate/(1024^2))
    @time run_simulation(sim, app, state)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate/(1024^2))
  end
end
