#======================================================================
script_name=$0                   # nome dello script
exec_file=$1                     # nome algoritmo di fattorizzazione
directory="."                    # directory dove è presente lo script
host_file=""                     # host file di open mpi
run="mpirun"                         # comando di esecuzione: es mpirun
options=""  # opzioni per il comando di esecuzione
#======================================================================

# vai alla cartella di lavoro
cd $directory

# function usage{

# Formato di output
# NUM1, T1, T2, ... TM
# NUM2, T1, T2, ... TM
# ...   ..  ..  ... ..
# NUMN, T1, T2, ... TM

# Fattorizza i numeri da $2 a $4 con passo $3
# Ogni fattorizzazione viene iterata $5 volte

#     nome_script  algorit_da_exec  da  passo  a     iteraz
# es: script.sh    fattor           1   2      1000  50
#   N = (1000 - 1)/2
#   M = 50 

# }

for N in $(seq $2 $3 $4)
do
    # "NUM1 "
    echo -n "$N"  
    for M in $(seq 1 1 $5)
    do

	echo -n ", "		
	$run $options $exec_file $N

    done
	
    # "\n"
    echo -e "\n"
done
