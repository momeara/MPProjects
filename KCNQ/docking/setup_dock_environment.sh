#!/bin/bash



# supported cluster types:
# LOCAL, SGE, SLURM
export CLUSTER_TYPE=SLURM

export SLURM_ACCOUNT=maom0
export SLURM_MAIL_USER=maom@umich.edu
export SLURM_MAIL_TYPE=BEGIN,END
export SLURM_PARTITION=standard

export SCRATCH_DIR=/scratch/maom_root/maom0/maom

export DOCK_TEMPLATE=/home/maom/opt/dock_campaign_template

export DOCKBASE=/home/maom/opt/DOCK37
export PATH="${PATH}:${DOCKBASE}/bin:${DOCK_TEMPLATE}/scripts"

export OMEGA_ENERGY_WINDOW=12
export OMEGA_MAX_CONFS=600


conda activate dock

