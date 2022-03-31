using SolitonBEC
using Printf
## Simulation
std_sim = Simulation(20e-9,
  1e-12,
  50000,
  100)
## Apparata
std_apparatus = Apparatus(1,
               1,
               1, 
               1, 
               1)
## InitialStates      
InitialState1 =  InitialState(0.001,
                20e-12,
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
