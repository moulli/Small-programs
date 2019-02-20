close all; %clc



%% Trying to compute own sparse matrix for saving:

ztemp = zgridGeo01;
[lx, ly, lz, ld] = size(ztemp.Zcorrelations);

findind = find(ztemp.Zcorrelations ~= 0);
lfind = length(findind);
condense = cell(lfind, 4);
for i = 1:lfind
    [ix, iy, iz, id] = ind2sub(size(ztemp.Zcorrelations), findind(i));
    condense{i, 1} = findind(i);
    condense{i, 2} = ztemp.Zneurons{ix, iy, iz, id};
    condense{i, 3} = ztemp.Zcorrelations(ix, iy, iz, id);
    condense{i, 4} = ztemp.Zneuron_number(ix, iy, iz, id);
    if mod(i, floor(lfind/20)) == 0
        fprintf('%.0f %% done. \n', 100*i/lfind);
    end
end

scond = struct;
castok = 0;
if castok == 1
    scond.size = cast(size(ztemp.Zcorrelations), 'single');
    scond.index = cast(findind, 'single');
    scond.correlations = cast(ztemp.Zcorrelations(findind), 'single');
    scond.number = cast(ztemp.Zneuron_number(findind), 'single');
else
    scond.size = size(ztemp.Zcorrelations);
    scond.index = findind;
    scond.correlations = ztemp.Zcorrelations(findind);
    scond.number = ztemp.Zneuron_number(findind);
end
maxcel = 0;
laycell = ztemp.Zneurons(findind);
for i = 1:lfind
    maxcel = max([maxcel, length(laycell{i})]);
end
neu_temp = zeros(lfind, maxcel);
for i = 1:lfind
    neu_temp(i, 1:length(laycell{i})) = laycell{i}';
end
scond.neurons = sparse(neu_temp);
ZTEMP = struct;
ZTEMP.method = ztemp.method;
ZTEMP.names = ztemp.names;
ZTEMP.paths = ztemp.paths;
ZTEMP.comments = ztemp.comments;
ZTEMP.increment = ztemp.increment;
ZTEMP.xgrid = ztemp.xgrid;
ZTEMP.ygrid = ztemp.ygrid;
ZTEMP.zgrid = ztemp.zgrid;
ZTEMP.corvect = ztemp.Zcorvect;
ZTEMP.stru = scond;


% condense = cell(sum(ztemp.Zcorrelations ~= 0), 2);
% cindex = 1;
% for id = 1:ld
%     for ix = 1:lx
%         for iy = 1:ly
%             for iz = 1:lz
%                 if ~isempty(ztemp.Zcorrelations(ix, iy, iz, id))
%                     condense{cindex, 1} = [ix, iy, iz, id];
%                     condense{cindex, 2} = {ztemp.Zneurons{ix, iy, iz, id}, ...
%                         ztemp.Zcorrelations(ix, iy, iz, id), ztemp.Zneuron_number(ix, iy, iz, id)};
%                     cindex = cindex + 1;
%                     fprintf('
%                 end
%             end
%         end
%     end
% end
