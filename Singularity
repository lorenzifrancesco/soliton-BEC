Bootstrap: docker
From: julia:1.7.2

%post
    apt update -y
    apt upgrade -y
    apt install git -y
        
%runscript
    julia --threads=auto install.jl