Bootstrap:docker
From:debian:latest

%labels
    MAINTAINER Peter Kruczkiewicz
    DESCRIPTION Singularity image containing all requirements for the peterk87/nf-iav-illumina pipeline
    VERSION 1.1.0

%environment
    export PATH=/opt/conda/envs/nf-iav-illumina-1.1.0/bin:/opt/conda/bin:$PATH
    export LANG=C.UTF-8
    export LC_ALL=C.UTF-8

%files
    environment.yml /

%post
    apt-get update --fix-missing
    apt-get install -y wget bzip2 ca-certificates curl git procps
    apt-get clean
    rm -rf /var/lib/apt/lists/*
    wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.10-Linux-x86_64.sh -O /miniconda.sh 
    /bin/bash /miniconda.sh -b -p /opt/conda
    rm /miniconda.sh
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
    export PATH=/opt/conda/bin:$PATH
    conda install conda=4.7.12
    conda env create -f /environment.yml
    conda clean -tipsy
