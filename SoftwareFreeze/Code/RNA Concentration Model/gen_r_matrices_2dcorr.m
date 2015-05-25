%valArray = [.01 .05 .1 .5 1 5 10 50 100 500];
valArray = [.01 .1 1 10 100];
n = numel(valArray);
rValMat = zeros(n,n,n);
expt = 0;

for i = 1:n %lambda
    for j = 1:n %gamma
        for k = 1:n %mu
            %for m = 1:n %delta
                % assign initial lam, gam, mu values
                lam0 = valArray(i);
                gam0 = valArray(j);
                mu0 = valArray(k);
                abundMat(i,j,k) = mu0*lam0/(lam0+gam0);
                %del0 = valArray(m);
                % now assume mu has volume dependence
                vArray = 0.1:.1:5;%sort(5*rand(50,1));
                expt = 0;
                for p = 1:numel(vArray)
                    vol = vArray(p);
                    %cell = getpdf(lam0*vol,gam0,mu0,del0);
                    cell = getpdf(lam0,gam0,mu0*vol,1);
                    
                    %[a cutoff] = max(diff(cell));
                    %[a b] = max(cell(cutoff:end));
                    
                    %rnaArray(1,p) = b+cutoff;
                    %%[a rna] = max(cell(10:end));
                    %%rnaArray(1,p) = rna;
                    
                    if expt == 0
                        expt = cell;
                    else
                        expt = appendCell(expt,cell);
                    end
                end
                %corr2(vArray,rnaArray)
                %scatter(vArray,rnaArray)
                %calculate correlation coefficient for expt
                % r =
                rValMat(i,j,k) = calcCorr(expt);
            %end
        end
        
    end
end

%plotting
[x1 y1 z1] = meshgrid(1:5,1:5,1:5);
scatter3(x1(:),y1(:),z1(:),50,rValMat(:),'filled');
xlabel('gamma');
ylabel('lambda');
zlabel('mu');

[x y] = meshgrid(1:5,1:5);
twoDArray = rValMat(1,:,:);
size(twoDArray)
f = squeeze(twoDArray)
figure;scatter(x(:),y(:),50,f(:),'filled')

%
load('mu_variable.mat');
muMat = rValMat;
load('lam_variable.mat');
lamMat = rValMat;
load('gam_variable.mat');
gamMat = rValMat;
load('del_variable.mat');
delMat = rValMat;