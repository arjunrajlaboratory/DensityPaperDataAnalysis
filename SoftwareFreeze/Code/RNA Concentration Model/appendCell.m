function [ expt ] = appendCell( expt_prev, new_cell )
% Append one pdf generated by getpdf  to a series of others

if size(new_cell,2) > size(expt_prev,2) % expand size of expt, add new cell
    expt = zeros(size(expt_prev,1),size(new_cell,2));
    expt(1:size(expt_prev,1),1:size(expt_prev,2)) = expt_prev;
    expt(end+1,:) = new_cell;
elseif size(new_cell,2) < size(expt_prev,2)
    tmp = zeros(1,size(expt_prev,2));
    tmp(1:size(new_cell,2)) = new_cell;
    expt = expt_prev;
    expt(end+1,:) = tmp;
else
    expt = expt_prev;
    expt(end+1,:) = new_cell;
end

end

