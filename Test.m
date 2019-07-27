load('Data369');
Geno=Geno369;
kt=2;   %Coefficient of k (number of crosses)
p=0.7;  %Scaling factor in the CBS method
n=100;  %Number of Simulation Replications
for i=1:n
    [Phenotype,values]=CBS(Geno,kt,p, eft, RF);
    resultPheno{1,i}=Phenotype;
    resultValue{1,i}=values;
end
RESULTPhenoCBS=cell2mat(resultPheno);
save RESULTPhenoCBS
RESULTValueCBS=cell2mat(resultValue);
save RESULTValueCBS
