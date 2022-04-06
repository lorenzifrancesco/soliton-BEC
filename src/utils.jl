function graphics_config()
    Plots.plotly()
    Plots.default(size=(1300, 800),
                  guidefont=("times", 10), 
                  tickfont=("times", 10))
  end

function mem_estimate(sim::Simulation)
    return 2* sim.T/sim.dt * sim.S/sim.ds * 128
end

function describe_simulation(sim::Simulation, app::Apparatus, state::InitialState)
    
end