function [tata tata_is_good bre bre_is_good inr inr_is_good dpe dpe_is_good] = promoterelements(pre,trans);

%sequence of each promoter element
%start location of each promoter element (from start of coding region = +1)

%pre = 'ccgcccaaccggcgtccgcctataaaaagctgagtgttgacgtcagcgtt';
%trans = 'CTCTTCCGCCGTCGTCGCCGCCATCCTCGGCGCGACTCGCTTCTTTCGGT';
trans = lower(trans);

%This first part identifies *potential* promoter elements, possibly without
%perfect sequence matches.

%TATA box
tlist = find(pre == 't');
tlist(find(tlist>numel(pre)-3))=[];
keep = false(size(tlist));
for i = 1:numel(tlist)
    %disp(pre(tlist(i):tlist(i)+3))
    if pre(tlist(i):tlist(i)+3) == 'tata'
        keep(i) = true;
    end
end
tlist = tlist(keep);
if numel(tlist) > 1
    [val,pos] = min(abs(tlist-19));
    tata_start = tlist(pos);
elseif numel(tlist) == 1
    tata_start = tlist;
else
    tata_start = 0;
    tata = [];
end

if tata_start ~= 0
    tata_end = tata_start + 3;
    while tata_end ~= numel(pre)
        if pre(tata_end+1) ~= 'a'
            break
        else
            tata_end = tata_end + 1;
        end
    end
tata = pre(tata_start:tata_end);
tata_start = tata_start-50;
end

%BRE
if tata_start > 7 & pre((tata_start-4):(tata_start-1))=='cgcc'
    bre = pre((tata_start-7):(tata_start-1));
    bre_start = tata_start-7;
else
    clist = find(pre == 'c');
    clist(find(clist<4))=[];
    keep = false(size(clist));
    for i = 1:numel(clist)
        if clist(i)+3<numel(pre) & (pre(clist(i):clist(i)+3) == 'cgcc')
            keep(i) = true;
        end
    end
    clist = clist(keep);
    if numel(clist) > 1
        [val,pos] = min(abs(clist-13));
        bre_start = clist(pos)-3;
        bre = pre(bre_start:bre_start+6);
    elseif numel(clist) == 1
        bre_start = clist-3;
        bre = pre(bre_start:bre_start+6);
    else
        bre_start = [];
        bre = [];
    end
end

%Inr
inr = [pre((end-1):end) trans(1:5)];
inr_start = 49;

%DPE
dpe_start = 28;
dpe = trans(28:32);

%***QC step***%

%TATA
if tata_start ~= 0
    tata_is_good = true;
else
    tata_is_good = false;
end

%BRE
if bre_start~=0 & find(bre(1) == ['g' 'c']) & find(bre(2) == ['g' 'c']) & ...
        find(bre(3) == ['g' 'a']);
    bre_is_good = true;
else
    bre_is_good = false;
end

%Inr
if inr_start~=0 & (inr(3) == 'a') & find(inr(1) == ['t' 'c']) & find(inr(2) == ['t' 'c']) ...
        & find(inr(5) == ['t' 'a']) & find(inr(6) == ['t' 'c']) & ...
        find(inr(7) == ['t' 'c'])
    inr_is_good = true;
else
    inr_is_good = false;
end

%DPE
if dpe_start~=0 & (dpe(2) == 'g') & find(dpe(1) == ['a' 'g']) & find(dpe(3) == ['a' 't']) ...
        find(dpe(4) == ['c' 't']) & find(dpe(5) == ['a' 'g' 'c'])
    dpe_is_good = true;
else
    dpe_is_good = false;
end

tata;
tata_start;
tata_is_good;

bre;
bre_start;
bre_is_good;

inr;
inr_start;
inr_is_good;

dpe;
dpe_start;
dpe_is_good;