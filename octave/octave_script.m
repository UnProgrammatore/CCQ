M=csvread('esempio');

X=M(:,1);
Y=[];
Error=[];

for n=1:size(M)(1)
   Error=[Error, std(M(n,2:size(M)(2)))];
   Y=[Y, mean(M(n,2:size(M)(2)))];
end

xlabel(' Factorized ');
ylabel(' Time ');
errorbar(Error,X,Y);


print('-dbmp', 'Costo Computazionale');
