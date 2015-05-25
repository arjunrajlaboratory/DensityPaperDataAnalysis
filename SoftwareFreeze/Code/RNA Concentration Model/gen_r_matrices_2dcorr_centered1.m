%valArray = [.01 .05 .1 .5 1 5 10 50 100 500];
valArray = [.01 .05 .1 .1 1 5 10 50 100 500];
n = numel(valArray);
rValMat = zeros(n,n,n);
expt = 0;

for i = 1:n %lambda
    text = sprintf('on i=%d of %d',i,n);
    disp(text);
    for j = 1:n %gamma
        for k = 1:n %mu
            %for m = 1:n %delta
                % assign initial lam, gam, mu values
                lam0 = valArray(i);
                gam0 = valArray(j);
                mu0 = valArray(k);
                %abundMat(i,j,k) = mu0*lam0/(lam0+gam0);
                %del0 = valArray(m);
                % now assume mu has volume dependence
                vArray = 0.3:.03:1.7;%sort(5*rand(50,1));
                expt = 0;
                for p = 1:numel(vArray)
                    vol = vArray(p);
                    %cell = getpdf(lam0*vol,gam0,mu0,del0);
                    cell = getpdf(lam0,gam0/vol,mu0,1);
                    
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
                [corr avgN] = calcCorr(expt);
                rValMat(i,j,k) = corr;
                abundMat(i,j,k) = avgN;
            %end
        end
        
    end
end

%plotting
[x1 y1 z1] = meshgrid(1:n,1:n,1:n);
scatter3(x1(:),y1(:),z1(:),50,abundMat(:),'filled');
xlabel('gamma');
ylabel('lambda');
zlabel('mu');
set(gca,'XTick',1:10)
set(gca,'YTick',1:10)
set(gca,'ZTick',1:10)
set(gca,'XTickLabel',{.01, .05, .1, .1, 1, 5, 10, 50, 100, 500})
set(gca,'YTickLabel',{.01, .05, .1, .1, 1, 5, 10, 50, 100, 500})
set(gca,'ZTickLabel',{.01, .05, .1, .1, 1, 5, 10, 50, 100, 500})
% 
% [x y] = meshgrid(1:5,1:5);
% twoDArray = rValMat(1,:,:);
% size(twoDArray)
% f = squeeze(twoDArray)
% figure;scatter(x(:),y(:),50,f(:),'filled')
% 
% %
% load('mu_variable.mat');
% muMat = rValMat;
% load('lam_variable.mat');
% lamMat = rValMat;
% load('gam_variable.mat');
% gamMat = rValMat;
% load('del_variable.mat');
% delMat = rValMat;