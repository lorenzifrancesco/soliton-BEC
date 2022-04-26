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
  50e-3, #T
  5e-5, #dt
  L, #S
  L * 1e-3, #ds
)


# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "NPSE",
  "ellipse",
)

#--------- Potential ---------
function elliptical(ecc::Float64)
  r = Potential(
    "ellipse", # type 
    5e-6, #width
    15e-6, #position
    5e-6, #energy
    ecc, #ϵ
    L / E(2 * pi, ecc) #a
  )
  return r
end

# --------- Apparata ---------
std_apparatus = Apparatus(
  6.941 * 1.660539e-27, #m (conversion AMU -> kg)
  0, # as
  2 * pi * 710, # ω_perp
  4e3, #N
  1.77e-11, #γ
  2 * pi * 4,# ω_z
)

# --------- InitialStates ---------   
InitialState1 = InitialState(
  "sech", #type
  1.7e-6, # width
  0e6 #v0
)

## Configurations
configs = []
for ecc in [0, 0.5, 0.75, 0.9]
  potential = elliptical(ecc)
  push!(configs, (num, khaykovich_gpe, potential, std_apparatus, InitialState1))
end

display(configs)
mem_limit = 15000000000 #byte
cnt = 1

pyplot()
plt_width = 800
plt_height = 600
## Curvature potential--------------------------------------------

fig1 = plot(title="",
  xlabel="s/L",
  ylabel="U_Q/E_a",
  size=(plt_width, plt_height))

for (num, sim, pot, app, state) in configs
  characteristic_energy = hbar^2 / (app.m * pot.a^2)

  space_steps = Int(floor(num.S / num.ds))
  space = LinRange(-num.S / 2, num.S / 2, space_steps)
  plot!(fig1, space / L, -abs.(map(s -> potential(sim, app, pot, s), space)) / characteristic_energy, reuse=false,
    label="ϵ = $(pot.ϵ)")
end


## Non Interacting ground state--------------------------------------------

fig2 = plot(title="",
  xlabel="s/L",
  ylabel="ρ(s) * a",
  reuse=false,
  size=(plt_width, plt_height))

for (num, sim, pot, app, state) in configs
  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit

    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))

    # Natural orientation choice for plot 
    coeffs = get_coefficients(sim, app, pot, state)

    time_steps = Int(floor(num.T / num.dt))
    space_steps = Int(floor(num.S / num.ds))
    time = LinRange(0, num.T, time_steps)
    space = LinRange(-num.S / 2, num.S / 2, space_steps)

    time, space, ψ = @time ground_state_solve(num, coeffs)
    plot!(fig2, space / L,
      abs.(ψ) .^ 2 * pot.a,
      label="ϵ = $(pot.ϵ)")

    display(fig1)
    display(fig2)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end

## Spontaneous Symmetry Breaking--------------------------------------------
configs_gamma = []

fixed_potential = elliptical(0.9)
mass = 6.941 * 1.660539e-27
ω_perp = 2 * pi * 710
Num = 4e3

as_apparatus(as) = Apparatus(
  mass, #m (conversion AMU -> kg)
  as, # as
  ω_perp, # ω_perp
  Num, #N
  1.77e-11, #γ
  2 * pi * 4,# ω_z
)
l_perp = sqrt(hbar / (mass * ω_perp))

for gamma in [-1, -1.5, -2, -3, 2, 3]
  as = gamma * l_perp^2 / (2 * fixed_potential.a * (Num - 1))
  display(as_apparatus(as))
  push!(configs_gamma, (num, khaykovich_gpe, fixed_potential, as_apparatus(as), InitialState1))
end
p = Plots.palette(:rainbow_bgyr_35_85_c72_n256, length(configs_gamma) + 3)


fig3 = plot(title="",
  xlabel="s/L",
  ylabel="ρ(s) * a",
  reuse=false,
  size=(plt_width, plt_height),
  palette=p)

for (num, sim, pot, app, state) in configs_gamma
  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit

    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))

    # Natural orientation choice for plot 
    coeffs = get_coefficients(sim, app, pot, state)

    time_steps = Int(floor(num.T / num.dt))
    space_steps = Int(floor(num.S / num.ds))
    time = LinRange(0, num.T, time_steps)
    space = LinRange(-num.S / 2, num.S / 2, space_steps)

    time, space, ψ = @time ground_state_solve(num, coeffs)
    γ = (2 * pot.a * app.as * (app.N - 1)) / (l_perp^2)
    plot!(fig3, space / L,
      abs.(ψ) .^ 2 * pot.a,
      label="γ = $(γ)")

    display(fig3)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end