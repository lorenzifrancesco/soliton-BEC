Bootstrap: docker
From: julia:1.7.2
%post
    apt update
    apt install git
    git clone https://github.com/lorenzifrancesco/soliton-BEC
    cd soliton-BEC
    git checkout transmission-exploration
        
%runscript
    julia -c ] activate SolitonBEC
