Bootstrap: docker
From: julia:1.7.2

%post
    apt update -y
    apt upgrade -y
    apt install git -y
    git clone https://github.com/lorenzifrancesco/soliton-BEC
    cd soliton-BEC
    git checkout 3D-GPE
    
%runscript
    julia --threads=auto install.jl