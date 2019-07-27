function [Phenotype,values] = CBS(Geno, kt, p,eft,RF)
rp = randperm(369);
Geno = Geno(:,:,rp(1:200)); %Select 200 indiviuals
HT = 10;  %Total number of generations
cprog = 10; %Cost of making a progeny
ccross = 5; %Cost of making a cross
budget = 20500*ones(HT+1,1); %Budget
AllPhenos=cell(HT+1,1);
AllPhenos{1} =eft' * reshape(sum(Geno,2),size(Geno,1),[]);
M3Phenos = [min(AllPhenos{1}), mean(AllPhenos{1}), max(AllPhenos{1})];
potentials = [2*min(min(Geno,[],2),[],3)'*eft, 2*max(max(Geno,[],2),[],3)'*eft];

figpheno = figure;
ax1 = axes('Parent',figpheno,'FontName','Times New Roman','FontSize',28);
hold(ax1,'on');
plotpheno(figpheno, AllPhenos, M3Phenos, potentials, 0);


for t = 1:HT
    selection = GSelect(Geno, RF, eft, AllPhenos{t}, cprog, ccross, budget(t), HT+1-t, kt, p)
    budget(t+1) = budget(t) - size(selection,1)*ccross - sum(selection(:,3)) * cprog;
    [Geno, AllPhenos{t+1}, M3Phenos(t+1,:), potentials(t+1,:)] = reproduce(Geno, RF, eft, selection);
    plotpheno(figpheno, AllPhenos, M3Phenos, potentials, t);
end
values=[potentials(1:t+1,2)'; M3Phenos(:,3)'; M3Phenos(:,2)'; M3Phenos(:,1)'; potentials(1:t+1,1)'; budget']
Phenotype=AllPhenos{end,1};
end

%Reproduction
function [G, pheno, m3pheno, potential] = reproduce(G, RF, eft, selection)
nc = size(selection,1);
BP = G(:,:,reshape(selection(:,1:2)',1,[]));
G = repmat(G(:,:,1),1,1,sum(selection(:,3)));
m = 0;
for c = 1:nc
    np = selection(c,3);
    G(:,:,m+[1:np]) = cross(BP(:,:,2*c-1),BP(:,:,2*c),RF,np);
    m = m + np;
end
pheno = eft' * reshape(sum(G,2),size(G,1),[]);
m3pheno = [min(pheno), mean(pheno), max(pheno)];
potential = [2*min(min(G,[],2),[],3)'*eft, 2*max(max(G,[],2),[],3)'*eft];
end
function Y = cross(L1,L2,RF,k)
RN = rand(size(L1,1),2,k);
RC = RN<=repmat([0.5; RF],1,2,k);
SD = (1-cumprod(1-RC*2))/2;
Y = repmat([L1(:,1),L2(:,1)],1,1,k);
PR = repmat([L1(:,2),L2(:,2)],1,1,k);
Y(SD==1) = PR(SD==1);
end

%Selection
function selection = GSelect(Geno, RF, eft, pheno, cprog, ccross, budget, T, kt, p)
k=(kt*T)-1;
[~, bestmale] = sort(pheno,'descend');
x=size(Geno,3);
u=zeros(1,x);
bestfemale=zeros(1,k);
for j=1:k
    for i=1:x
        u(i)=ComplementarityScore(Geno(:,:,bestmale(j)),Geno(:,:,i),eft,p);
    end
    [~,ind]=max(u);
    bestfemale(j)=ind;
end
np=20;
selection = [bestmale(1:k)', bestfemale(1:k)', np*ones(k,1)];
end
function [score]=ComplementarityScore(ind1,ind2,eft,p)
% p is a parameter between 0 and 1
add=sum(ind1,2)+sum(ind2,2);
M=(add.^p);
score=(eft'*M);
end

%Plotting
function plotpheno(figpheno, AllPhenos, M3Phenos, potentials, t)
figure(figpheno);
fill([0:t, t:-1:0], [M3Phenos(:,1)', M3Phenos(end:-1:1,3)'], [0.8, 0.8, 1],'edgecolor','none');
hold on;
for i = 1:t+1
    [nn, yy] = hist(AllPhenos{i});
    dy = mean(diff(yy))/2;
    xf = (i-1) + [reshape([zeros(1,10);nn;nn],[],1);0;0]/sum(nn);
    yf = [reshape([yy-dy;yy-dy;yy+dy],[],1);yy(10)+dy;yy(1)-dy];
    fill(xf,yf,[0 0 1]);
end
plot([0:t], M3Phenos(:,2)','-r');
fill([0:t,t,0]',potentials([1:t+1,1,1],1),0.4 + [0 0 0]);
fill([0:t,t,0]',potentials([1:t+1,1,1],2),0.4 + [0 0 0]);
grid on;
xlim([0 t+1]);
ylim([potentials(1,1) potentials(1,2)]);
xlabel('Generation number','FontName','Times New Roman','FontSize',28);
ylabel('Phenotype value','FontName','Times New Roman','FontSize',28);
drawnow
end
