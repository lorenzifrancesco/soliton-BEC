using SolitonBEC
using Printf

## Simulation
std_sim = Simulation(1e-3,
  1e-5,
  100e-6,
  5e-7)

## Apparata
std_apparatus = Apparatus(1e-26,
               1e-6,
               4e14, 
               1000, 
               5e-6)

## InitialStates      
InitialState1 =  InitialState(1,
                2e-12,
                0)    

## Configuration
configs = [(std_sim, std_apparatus, InitialState1),]
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
