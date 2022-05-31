using NLSolve
using Plots


function f(μ , N)
  return (1-  μ)^(3/2) - 3/2*(1-μ)^(1/2) +3/2/sqrt(2) * -0.21e-9 *  N / 1.4321331901586098e-6 = 0
end

initial_N = 4e3
n = LinRange(0, initial_N, 20)

solved_μ = nlsolve(x->f(x, n) , [0.5])

fig = plot(n, solved_μ)
display(fig)