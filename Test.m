load('Data369');
Geno=Geno369;
kt=2;
p=0.7;
n=10;
for i=1:n
    CBS(Geno,kt,p, eft, RF);
    resultPheno{1,i}=Phenotype;
    resultValue{1,i}=values;
end
RESULTPhenoCBS=cell2mat(resultPheno);
save RESULTPhenoCBS
RESULTValueCBS=cell2mat(resultValue);
save RESULTValueCBS