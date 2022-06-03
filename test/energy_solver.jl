using NLsolve
using Plots


function f!(F, μ , N)
  F .=  (1 .-  μ).^(3/2) .- 3/2*(1 .- μ).^(1/2) .- 3/2/sqrt(2) * 0.21e-9 *  N / 1.4321331901586098e-6
  return nothing
end

function j!(J, μ)
  J .= -3/2*(1 .-  μ).^(1/2) .- 3/4*(1 .- μ).^(-1/2)
  return nothing
end 

function g!(G, μ , gamma)
  G .= abs.(1 .-  μ).^(3/2) .- 3/2*abs.(1 .- μ).^(1/2) .+ 3/2/sqrt(2)*gamma
  return nothing
end

function gamma!(μ)
  return ((1 .-  μ).^(3/2) .- 3/2*(1 .- μ).^(1/2)) * 2/3*sqrt(2)
end

initial_N = 4e3
n_list = LinRange(0, initial_N, 20)
solved_μ = Array{Array{Float64}}(undef, length(n_list))
cnt = 1
print("\n Gamma values: \n")
print( 3/2/sqrt(2) * 0.21e-9 * n_list / 1.4321331901586098e-6)
fig = scatter()

gamma_list = LinRange(0, 1.20, 500)

solution_μ = Array{Array{Float64}}(undef, length(gamma_list))

for gamma in gamma_list
  global cnt
  g_func!(G, μ) = g!(G, μ, gamma)
  solution_μ[cnt] = nlsolve(g_func! ,  [0.0; 1.0]).zero
  scatter!([gamma], [nlsolve(g_func!,  [0.0; 1.0]).zero[1]])
  scatter!([gamma], [nlsolve(g_func!, [0.0; 1.0]).zero[2]])
  cnt += 1
end
display(fig)
# for gamma in gamma_list
#   global cnt
#   target_error = 0.001
#   lower = -1.0
#   upper = 1.0
#   let error
#       let midpoint
#           error = target_error + 1
#           while error > target_error
#               midpoint = (lower + upper) / 2
#               if f(midpoint, gamma) > 0
#                   lower = midpoint
#               else
#                   upper = midpoint
#               end
#               print(midpoint)
#               error = abs(f(midpoint, gamma))
#           end
#           solution_μ[cnt] = midpoint
#       end
#   end
#   cnt+=1
# end
print("\nSolution\n")
print(solution_μ)
display(fig)

μ_values = LinRange(0, 1, 100)
# fig2 = plot(μ_values, gamma(μ_values))
# display(fig2)