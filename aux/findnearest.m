function idx = findnearest(target,list)
% FINDNEAREST returns index of element in list closest to target
list = abs(list-target);
[~,idx] = min(list);
end