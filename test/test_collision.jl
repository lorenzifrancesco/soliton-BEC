using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)

# --------- Numerics ---------
L = 200e-6

num = Numerics(
  100e-4, #T
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
    L * 1e-3 / 3, #width
    40e-6, #position
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
  2.4e-6, # width
  velocity #v0
)

## Configurations
configs = []
for energy in LinRange(1e-37, 30e-37, 5)
  potential = barrier_height(energy)
  push!(configs, (num, khaykovich_gpe, std_apparatus, potential, InitialState1))
end

display(configs)
mem_limit = 15000000000 #byte
cnt = 1

pyplot()
p = Plots.palette(:rainbow_bgyr_35_85_c72_n256, length(configs) + 3)

## Soliton - barrier collision --------------------------------------------
@time run_dynamics(configs[1]...)

fig1 = plot(title="|ψ|^2 after collision with barrier",
  xlabel="space [mm]",
  ylabel="|ψ|^2",
  reuse=false,
  size=(800, 400))

for (num, sim, app, pot, state) in configs

  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))

    coeffs = get_coefficients(sim, app, pot, state)

    time, space, ψ, ψ_spect = @time ssfm_solve(num, coeffs)
    plot!(space * 1e3,
      abs.(ψ[:, end]) .^ 2,
      label="peak potential = $(pot.energy / (sqrt(2 * pi) * pot.width) /1.602176634e-19)   [eV]",
      lw=1,
      color=p[cnt])

    plot!(space * 1e3,
      abs.(map(s -> potential(sim, app, pot, s), space)) / hbar,
      lw=1,
      ls=:dot,
      label=nothing,
      color=p[cnt])
    display(space)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end
display(fig1)
