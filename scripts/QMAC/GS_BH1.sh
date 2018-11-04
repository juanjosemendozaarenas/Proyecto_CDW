#!/bin/bash
# Submits a single job to the QMAC cluster
# The name of the input and output files are given by the argument to this script

# Combine stdout and sterr into one output file:
#$ -j y
# Use "bash" shell:
#$ -S /bin/bash
# Change the output directory for stdout:
#$ -o sge_output
# Name the job:
#$ -N GS_BH
# Use current directory as working root:
#$ -cwd
# Send job notifications by email
#$ -m beas
# Email address to send notifications to (change to your address)
#$ -M JuanJose.Mendoza-Arenas@physics.ox.ac.uk

# Make directory for output
mkdir --parents sge_output

# Set the environment variables required for the simulation
. /share/apps/tnt/scripts/set_tnt_vars.sh

# Run the application
../bin/Ground_State -d ../output/ -o ../initfiles/initial_GS_BH1.mat
