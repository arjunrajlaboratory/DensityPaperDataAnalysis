valArray = [.01 .05 .1 .1 1 5 10 50 100 500];
n = numel(valArray);
rValMat = zeros(n,n);
expt = 0;

for i = 1:n %p1
    text = sprintf('on i=%d of %d',i,n);
    disp(text);
    for j = 1:n %p2
                p1 = valArray(i);
                p2 = valArray(j);
                %abundMat(i,j) = mu0*lam0/(lam0+gam0);
                vArray = 0.3:.03:1.7;
                expt = 0;
                for p = 1:numel(vArray)
                    vol = vArray(p);
                    cell = calc_burst_pdf(p1*vol,p2);
                    
                    if expt == 0
                        expt = cell;
                    else
                        expt = appendCell(expt,cell);
                    end
                end
                rValMat(i,j) = calcCorr(expt);
    end
end

%plotting
[x1 y1] = meshgrid(1:n,1:n);
scatter(x1(:),y1(:),50,rValMat(:),'filled');
xlabel('lambda/delta');
ylabel('mu/gamma');
set(gca,'XTick',1:10)
set(gca,'YTick',1:10)
set(gca,'XTickLabel',{.01, .05, .1, .1, 1, 5, 10, 50, 100, 500})
set(gca,'YTickLabel',{.01, .05, .1, .1, 1, 5, 10, 50, 100, 500})
%xlabel('gamma');
%ylabel('lambda');
%zlabel('mu');
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