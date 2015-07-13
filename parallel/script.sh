#!/bin/bash

#PBS -A INFNG_test 
#PBS -l select=1:ncpus=1:mpiprocs=1+16:ncpus=16:mem=100gb:mpiprocs=1 
#PBS -l walltime=24:00:00 

#PBS -N CCQi-n16-50l

EXE=executables/qs-intel
module load intel intelmpi
#EXE=executables/qs-gnu
#module load gnu openmpi
#export GOMP_CPU_AFFINITY=0-15

########################################################
cd /galileo/home/userexternal/ralfieri/git/CCQ/parallel/
cat $PBS_NODEFILE
date
########################################################

# echo "#20 cifre (base_dim = 148)"   ; mpirun $EXE  1 18567078082619935259 4000 10000000 10 1000 50
# echo "#30 cifre (base_dim = 222)"   ; mpirun $EXE  1 350243405507562291174415825999 9000 10000000 10 1000 50
# echo "#40 cifre (base_dim = 880)"   ; mpirun $EXE  1 5705979550618670446308578858542675373983 40000 100000000 10 10000 50
# echo "#45 cifre (base_dim = 1633)"  ; mpirun $EXE  1 732197471686198597184965476425281169401188191 60000 100000000 10 10000 100
 echo "#50 cifre (base_dim = 2540)"  ; mpirun $EXE  1 53468946676763197941455249471721044636943883361749 120000 100000000 10 10000 100
# echo "#55 cifre (base_dim = 3422)"  ; mpirun $EXE  1 5945326581537513157038636316967257854322393895035230547 180000 100000000 10 10000 100
# echo "#60 cifre (base_dim = 4842)"  ; mpirun $EXE  1 676292275716558246502605230897191366469551764092181362779759 250000 100000000 10 10000 100  
# echo "#65 cifre (base_dim = 6980)"  ; mpirun $EXE  1 530332320518131653109331688115692193258351283222071577636073763 400000 100000000 10 10000 100 
# echo "#71 cifre (base_dim = 8983)"  ; mpirun $EXE  1 1630245464892823711538270142482769515544228277190451217281557216919689 700000 100000000 10 10000 100 
# echo "#80 cifre (base_dim = 20902)" ; mpirun $EXE  1 3521851118865011044136429217528930691441965435121409905222808922963363310303627 1000000 100000000 10 10000 100


## INPUT: [n1] [n2] [dim_crivello_eratostene] [dim_intervallo] [dim_blocco]
## n1, n2: a scopo di debug l'applicazione prende in input due numeri che moltiplica tra loro producendo il numero
##         N da fattorizzare. Porre uno dei due numeri a 1 per inserire direttamente l'N da fattorizzare.
##
## dim_crivello_eratostene: per calcolare la base di fattori si calcolano prima i numeri primi da 
##                          0 a dim_crivello_eratostene (nota: la base di fattori sarà molto più piccola
##			    di questo parametro)
## dim_intervallo: dimensione insieme di A per i queali ogni salve fattorizza i rispettivi Q(A)
##
## dim_blocco: dato un intervallo ad uno slave questo lo separa in altri sottointervalli che assegna ad
##             ogni thread. ogni thread calcola dim_blocco valori alla volta e spedisce al master i risultati.
## fact_print: numero che indica ogni quante fattorizzazioni ricevute dal master stampare un feedback in stdout.


