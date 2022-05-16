using SolitonBEC
using Printf
using Plots
using Elliptic
hbar = 6.62607015e-34 / (2 * pi)
width = 1.42e-6
mass = 6.941 * 1.660539e-27
omega_perp = 2 * pi * 710
l_perp = sqrt(hbar / (mass * omega_perp))
N = 4e3
as = -0.21e-9
interaction_g = 2 * hbar * omega_perp * abs(as)

# mass as total mass or mass per atom??
time_unit = hbar^2/(mass * interaction_g * N)
space_unit = hbar^3/(mass * interaction_g^2 * N^2)
velocity_unit = interaction_g * N / hbar

# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "GPE",
  "barrier",
)

# --------- Apparata ---------
std_apparatus = Apparatus(
  mass, #m (conversion AMU -> kg)
  as, # as
  omega_perp, # ω_perp
  N, #N
  1.77e-11, #γ
  2 * pi * 4,# ω_z
)

# Energy = 1 / 2 * (std_apparatus.as^2 / (width^2) + std_apparatus.as^2 / (2 * l_perp^2) + width^2 / std_apparatus.as^2 - interaction_g * std_apparatus.as^3 / (l_perp * width^2)) * hbar * omega_perp


# --------- Numerics ---------
T = 2 * time_unit
L = 2 * space_unit
print("\ntime_unit: ", time_unit)
print("\nspace_unit: ", space_unit)
num = Numerics(
  T, #T
  T * 1e-3, #dt
  L, #S
  L * 1e-3, #ds
)

v0 = 1.5e6
# --------- InitialStates ---------   
InitialState1 = InitialState(
  "sech", #type
  width, # width
  v0 #v0
)

#--------- Potential ---------

function barrier_height(energy::Float64)
  r = Potential(
    "barrier", # type 
    L * 1e-3 * 10, #width
    30e-6, #position
    energy, #energy
    0, #ϵ
    10e-6 #a
  )
  return r
end


## Configurations
configs = []
for energy in LinRange(1e-37, 100e-37, 2)
  potential = barrier_height(energy)
  push!(configs, (num, khaykovich_gpe, std_apparatus, potential, InitialState1))
end

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
  size=(800, 400),
  legend=:topleft)

for (num, sim, app, pot, state) in configs

  global cnt
  @printf("\n------Running simulation # %3i \n", cnt)
  estimate = mem_estimate(num)

  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate / (1024^2))

    coeffs = get_coefficients(sim, app, pot, state)

    time, space, ψ, ψ_spect = @time ssfm_solve(num, coeffs)
    # plot!(space * 1e3,
    #   abs.(ψ[:, end]) .^ 2,
    #   label="peak potential = $(pot.energy / (sqrt(2 * pi) * pot.width) /1.602176634e-19)   [eV]",
    #   lw=1,
    #   color=p[cnt])

    # plot!(space * 1e3,
    #   abs.(map(s -> potential(sim, app, pot, s), space)) / hbar,
    #   lw=1,
    #   ls=:dot,
    #   label=nothing,
    #   color=p[cnt])
    # display(space)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate / (1024^2))
  end

  cnt += 1

end