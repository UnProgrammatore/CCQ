#           da passo a    outfile_prime     outfile_eratos
# script.sh 1  2     1000 prim-1-2-1000.out erat-1-2-1000.out
#

echo "" > $1
echo "" > $2

#for N in $(seq $1 $2 $3)
#do
#    ./crivello_eratostene_ser $N >> $5
#    ./crivello_eratostene_par $N >> $4
#done

file=`cat sample.txt`

for N in $file
do
    ./crivello_eratostene_ser $N >> $2
    ./crivello_eratostene_par $N >> $1
done


