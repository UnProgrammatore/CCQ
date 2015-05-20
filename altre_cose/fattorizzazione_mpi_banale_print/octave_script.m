M=csvread('prova');

X=M(:,1);
Y=[];
Error=[];

for n=1:size(M)(1)
   Error=[Error, std(M(n,2:size(M)(2)))];
   Y=[Y, mean(M(n,2:size(M)(2)))];
end

errorbar(X,Y,Error);
xlabel(' Factorized ');
ylabel(' Time ');

print('-dbmp', 'Costo Computazionale');
