clear pdfplot;
lam = 1;
gam = 0.1;
mu = 50;
del = 4;

for i=1:100
    pdf = getpdf(lam,gam,mu,del/i);
    if i == 1 
        pdfplot = pdf;
    end
    
    if size(pdf,2) > size(pdfplot,2)
        tmp = pdfplot;
        pdfplot = zeros(size(pdfplot,1),size(pdf,2));
        pdfplot(1:size(tmp,1),1:size(tmp,2)) = tmp;
        pdfplot(i,:) = pdf;
    elseif size(pdf,2) < size(pdfplot,2)
        tmp = zeros(1,size(pdfplot,2))
        tmp(size(pdf)) = pdf;
        pdfplot(i,:) = tmp;
    else
        pdfplot(i,:) = pdf;
    end
end

imagesc(pdfplot)