### TEST FOR THE GROUND STATE OF ELLIPTICAL WAVEGUIDE FOR DIFFERENT ECCENTRICITY
## PLOTTING RESULTS IN NORMALIZED COORDINATES AS Salasnich - SciPost Physics (2022)

using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)
width = 1.42e-6
mass = 6.941 * 1.660539e-27
omega_perp = 2 * pi * 710
l_perp = sqrt(hbar / (mass * omega_perp))

# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "GPE",
  "barrier",
)

# --------- Apparata ---------
std_apparatus = Apparatus(
  mass, #m (conversion AMU -> kg)
  -0.21e-9, # as
  omega_perp, # ω_perp
  4e3, #N
  1.77e-11, #γ
  2 * pi * 4,# ω_z
)

interaction_g = 2 * hbar^2 / mass * std_apparatus.as * (std_apparatus.N - 1) / width^2 # 2/width^2
display(interaction_g)

# Energy = 1 / 2 * (std_apparatus.as^2 / (width^2) + std_apparatus.as^2 / (2 * l_perp^2) + width^2 / std_apparatus.as^2 - interaction_g * std_apparatus.as^3 / (l_perp * width^2)) * hbar * omega_perp

# adimensionalization
velocity = 0.8 * width
μ = interaction_g^2 / 8
Energy = (velocity^2 / (2 * mass) * std_apparatus.m + μ) * (std_apparatus.N)
display(Energy)
velocity_coefficient = Energy / (hbar * velocity)
display(velocity)
display(velocity_coefficient)

# --------- Numerics ---------
T = 10e-4
L = 2 * velocity * T

num = Numerics(
  T, #T
  T * 1e-3, #dt
  L, #S
  L * 1e-3, #ds
)

# --------- InitialStates ---------   
display(width)
InitialState1 = InitialState(
  "sech", #type
  width, # width
  velocity_coefficient #v0
)
#--------- Potential ---------

function barrier_height(energy::Float64)
  r = Potential(
    "barrier", # type 
    L * 1e-3 / 3, #width
    L * 3 / 4, #position
    energy, #energy
    0, #ϵ
    10e-6 #a
  )
  return r
end

## Configurations
configs = []

for energy_vs_g in LinRange(-0.5, -1.5, 4)
  energy = abs(interaction_g) * energy_vs_g
  potential = barrier_height(energy)
  push!(configs, (num, khaykovich_gpe, std_apparatus, potential, InitialState1))
end

#display(configs)
mem_limit = 15000000000 #byt

pyplot()
plt_width = 800
plt_height = 600
p = Plots.palette(:rainbow_bgyr_35_85_c72_n256, length(configs) + 3)

## Soliton - barrier collision --------------------------------------------
@time run_dynamics(configs[1]...)

cnt = 1
T = zeros(Float64, 1, length(configs))

# num stays constant
estimate = mem_estimate(num)
if estimate < mem_limit
  @printf("Estimated memory per simulation consumption: %4.1f MiB\n", estimate / (1024^2))
else
  @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
end

for (num, sim, app, pot, state) in configs
  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  # potential space index
  coeffs = get_coefficients(sim, app, pot, state)
  time, space, ψ, ψ_spect = @time ssfm_solve(num, coeffs)
  potential_idx = Int64(floor((pot.position+num.S/2) / num.ds)

  T[cnt] = sum(abs.(ψ[potential_idx:length(space)]) .^ 2 * num.ds)
  cnt += 1

end


# plot!(space * 1e3,
#   abs.(ψ[:, end]) .^ 2,
#   label="energy = $(pot.energy)",
#   lw=2,
#   color=p[cnt])

fig1 = plot(title="|ψ|^2 after collision with barrier",
  xlabel="space [mm]",
  ylabel="|ψ|^2",
  reuse=false,
  size=(plt_width, plt_height))
plot!(1:cnt,
  T,
  lw=1,
  ls=:dot,
  label=nothing,
  color=p[cnt])
display(fig1)
