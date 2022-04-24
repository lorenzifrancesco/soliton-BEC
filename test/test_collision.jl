### TEST FOR THE GROUND STATE OF ELLIPTICAL WAVEGUIDE FOR DIFFERENT ECCENTRICITY
## PLOTTING RESULTS IN NORMALIZED COORDINATES AS Salasnich - SciPost Physics (2022)

using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)

# --------- Numerics ---------
L = 100e-6

num = Numerics(
  50e-4, #T
  5e-6, #dt
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
    1e-7, #width
    15e-6, #position
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
velocity = 1.5e6
InitialState1 = InitialState(
  "sech", #type
  5e-6, # width
  velocity #v0
)

## Configurations
configs = []
for energy in LinRange(1e6, 40e6, 4)
  potential = barrier_height(energy)
  push!(configs, (num, khaykovich_gpe, std_apparatus, potential, InitialState1))
end

display(configs)
mem_limit = 15000000000 #byte
cnt = 1

pyplot()
p = Plots.palette(:jet1, 11)

## Soliton - barrier collision --------------------------------------------
#@time run_dynamics(configs[4]...)

fig1 = plot(title="|ψ|^2 after collision with barrier",
  xlabel="space [mm]",
  ylabel="|ψ|^2",
  reuse=false)

for (num, sim, app, pot, state) in configs

  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))

    coeffs = get_coefficients(sim, app, pot, state)

    time_steps = Int(floor(num.T / num.dt))
    space_steps = Int(floor(num.S / num.ds))
    time = LinRange(0, num.T, time_steps)
    space = LinRange(-num.S / 2, num.S / 2, space_steps)

    # plot(space / L, -abs.(map(s -> potential(sim, app, pot, s), space)), show=true, reuse=false, title="potential")
    space, time, ψ, ψ_spect = @time ssfm_solve(num, coeffs)
    plot!(fig1, space * 1e3,
      abs.(ψ[:, time_steps]) .^ 2,
      label="energy = $(pot.energy)",
      lw=2,
      palette=p)


  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end
display(fig1)
