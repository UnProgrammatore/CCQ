#cd /home/didattica/andrea.segalini/progetto

for N in $(seq $1 $2 $3)
do
	mpirun fattor $N
done

#mpirun --hostfile hosts hello
