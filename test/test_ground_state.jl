using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)

### SIMULATE GENERATION OF SOLITON IN Li-7 CONDENSATE

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

barrier = Potential(
  "barrier", #type
  5e-6, #width
  15e-6, #position
  5e-6, #energy
  0.9, #ϵ
  10e-6 #a
)

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
InitialState1 = InitialState(
  "sech", #type
  1.7e-6 # width
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


## Ground state--------------------------------------------

fig2 = plot(title="stationary |ψ|^2",
  xlabel="space [mm]",
  ylabel="ρ(s) * a",
  reuse=false)

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

    space, time, ψ = @time ground_state_solve(num, coeffs)
    plot!(fig2, space * 1e3,
      abs.(ψ) .^ 2 * pot.a,
      label="ϵ = $(pot.ϵ)")


    display(fig1)
    display(fig2)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end