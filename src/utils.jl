function config()
    Plots.plotly()
    Plots.default(size=(1300, 800),
                  guidefont=("times", 10), 
                  legendfont=("times", 12),
                  tickfont=("times", 10)
                  )
  end