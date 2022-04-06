using SolitonBEC
using Printf

### SIMULATE GENERATION OF SOLITON IN Li-7 CONDENSATE

# --------- Simulation ---------
khaykovich_gpe = Simulation(20e-3, #T
  5e-5, #dt
  1e-4, #S
  1e-7, #ds
  "GPE") #ds

khaykovich_npse = Simulation(10e-3, #T
  1e-5, #dt
  1e-4, #S
  1e-7, #ds
  "NPSE") #ds

npse_sim = Simulation(300e-3, #T      ##exibits Talbot
  5e-4, #dt
  100e-6, #S
  1e-8, #ds
  "NPSE") 
  
# --------- Apparata ---------
std_apparatus = Apparatus(6.941 * 1.660539e-27, #m (conversion AMU -> kg)
               -0.21e-9, #as
               2*pi*710, # ω_perp
               4e3, #N
               1.77e-11, # γ
               2*pi*4) # ω_z

## --------- InitialStates ---------   
InitialState1 =  InitialState(1.7e-6,
                0)    

## Configurations
configs = [(khaykovich_gpe, std_apparatus, InitialState1), 
         #(khaykovich_npse, std_apparatus, InitialState1), 
         ]
mem_limit = 15000000000 #byte
cnt=1
for (sim, app, state) in configs
  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(sim)
  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate/(1024^2))
    @time run_simulation(sim, app, state)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate/(1024^2))
  end
  cnt+=1
end