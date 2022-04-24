### TEST FOR THE GROUND STATE OF ELLIPTICAL WAVEGUIDE FOR DIFFERENT ECCENTRICITY
## PLOTTING RESULTS IN NORMALIZED COORDINATES AS Salasnich - SciPost Physics (2022)

using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)

# --------- Numerics ---------
L = 10e-6

num = Numerics(
  50e-3, #T
  5e-5, #dt
  L, #S
  L * 1e-3, #ds
)


# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "GPE",
  "barrier",
)

#--------- Potential ---------

function barrier_height(energy::Float64)
  r = Potential(
    "barrier", # type 
    L * 5e-2, #width
    L / 3, #position
    energy, #energy
    0, #ϵ
    10e-6 #a
  )
  return r
end


# --------- Apparata ---------
std_apparatus = Apparatus(
  6.941 * 1.660539e-27, #m (conversion AMU -> kg)
  -0.21e-9, # as
  2 * pi * 710, # ω_perp
  4e3, #N
  1.77e-11, #γ
  2 * pi * 4,# ω_z
)

# --------- InitialStates ---------   
velocity
InitialState1 = InitialState(
  "sech", #type
  1.7e-6, # width
  1.5e6 #v0
)

## Configurations
configs = []
for energy in 1:10
  energy = energy * 1e-6
  potential = barrier_height(energy)
  push!(configs, (num, khaykovich_gpe, potential, std_apparatus, InitialState1))
end

display(configs)
mem_limit = 15000000000 #byte
cnt = 1

pyplot()

## Curvature potential--------------------------------------------

fig1 = plot(title="curvature potential",
  xlabel="s/L",)

for (num, sim, pot, app, state) in configs
  characteristic_energy = hbar^2 / (app.m * pot.a^2)

  space_steps = Int(floor(num.S / num.ds))
  space = LinRange(-num.S / 2, num.S / 2, space_steps)
  plot!(fig1, space / L, -abs.(map(s -> potential(sim, app, pot, s), space)) / characteristic_energy, reuse=false,
    label="ϵ = $(pot.ϵ)")
end


## Soliton - barrier collision --------------------------------------------

for (num, sim, pot, app, state) in configs

  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))
    @time run_dynamics(num, sim, app, pot, state)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end