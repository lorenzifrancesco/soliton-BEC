attenuative = Apparatus(20e-9,
               1e-12,
               50000,
               100,
               0.2,
               0e-27,
               0.00e-3)
dispersive = Apparatus(2e-9,
               1e-13,
               77500,
               100,
               0.0,
               20e-27,
               0.00e-3)
nonlinear = Apparatus(20e-9,
               1e-12,
               27000,
               100,
               0.0,
               0e-27,
               3.00e-3)
solitonic = Apparatus(20e-9,
               1e-12,
               27000,
               100,
               0.0,
               -4e-27,
               3.00e-3)
smf_28 =Apparatus(10e-9,
              1e-13,
              55500,
              100,
              0.18,
              -21.58e-27,
              0.78e-3)
hnlf = Apparatus(1e-9,
            1e-13,
            5300,
            1,
            0.81,
            -7.09e-27,
            10.68e-3)
## InitialStates      
InitialState1 =  InitialState(0.001,
                20e-12,
                0)
InitialState2 =  InitialState(0.1,
                20e-12,
                0)
InitialState3 =  InitialState(0.46,
                7e-12,
                1)
InitialState4 =  InitialState(0.152,
                2.09e-12,
                1)  
high_power =  InitialState(0.152*9,
                2.09e-12,
                1)         
soliton = InitialState(3.4,
                20e-12,
                1)       
## Configuration
configs = [(smf_28, InitialState3), # first order soliton (attenuation limited)
          (hnlf, InitialState4)] # third order soliton (attenuation limited)

mem_limit = 15000000000 #byte

for (Apparatus, InitialState) in configs
  estimate = mem_estimate(Apparatus)
  if estimate < mem_limit
    @printf("Estimated memory consumption: %4.1f MiB\n", estimate/(1024^2))
    
    @time solve(Apparatus, InitialState)
  else
    @printf("Estimated memory consumption (%4.1f MiB) exceed maximum!\n", estimate/(1024^2))
  end
end
