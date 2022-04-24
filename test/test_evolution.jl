using SolitonBEC
using Printf
using Plots
using Elliptic
### SIMULATE GENERATION OF SOLITON IN Li-7 CONDENSATE

# --------- Numerics ---------
std = Numerics(
  T=50e-3,
  dt=5e-5,
  S=10e-5 * E(2 * pi, 0.9),
  ds=10e-8 * E(2 * pi, 0.9),
)

# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "GPE",
  "barrier",
) #ds

#--------- Potential ---------
ellipse = Potential(
  type="ellipse",)

barrier = Potential(
  type="barrier",)

# --------- Apparata ---------
std_apparatus = Apparatus(
  m=6.941 * 1.660539e-27, #m (conversion AMU -> kg)
  as=-0.21e-9, # m
  ω_perp=2 * pi * 710, # rad/s
  N=4e3,
  γ=1.77e-11,
  ω_z=2 * pi * 4,
  ϵ=0.9,
  a=10e-6
)

# --------- InitialStates ---------   
InitialState1 = InitialState(
  type="sech",
  width=1.7e-6
)

## Configurations
configs = [
  (khaykovich_gpe, std_apparatus, InitialState1),
  #(khaykovich_npse, std_apparatus, InitialState1), 
]
mem_limit = 15000000000 #byte
cnt = 1

angle = LinRange(0, 2 * pi, 100)

#plot(angle, -(sqrt(1 - 0.9^2) ./ (sin.(angle) .^ 2 + sqrt(1 - 0.9^2) * cos.(angle) .^ 2) .^ (3 / 2)) .^ 2, show=true)
for (sim, app, state) in configs
  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(sim)
  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))
    @time run_dynamics(num, sim, app, potential, state)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end
  cnt += 1
end