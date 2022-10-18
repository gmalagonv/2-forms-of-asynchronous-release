function all = normByRadius3(all)%,field)
%% function to normalize dits2cent by the following calculated ratio: 

if length([all.perSyn.synID])== length(unique([all.perSyn.synID])) && length([all.perSyn.synID]) == length(unique([all.allEv.synID]))
    disp('IDs OK')
else
    all = corrDuplisynIDS(all);
end

%% syn
if ~isfield(all.allEv ,'areaSynch')
    all = addFieldSyn2dEvs(all, 'areaSynch');
end


synIDs = unique([all.allEv.synID], 'stable');
numSyns = length(synIDs);
indxAll = (1:1:length(all.allEv));
for syn = 1 : numSyns
    
    radius = sqrt((all.perSyn(syn).areaSynch* 100000)/pi);
    mskEvs = ismember([all.allEv.synID], synIDs(syn));
    indxEvSyn = indxAll(mskEvs);
    mskSynch =  [all.allEv(indxEvSyn).synch];
    maxdist2centFSynch =  max([all.allEv(indxEvSyn(mskSynch)).dist2centFSynch]);
    
    for ev = 1:length(indxEvSyn)
        all.allEv(indxEvSyn(ev)).norm3 = all.allEv(indxEvSyn(ev)).dist2centFSynch/radius;
        all.allEv(indxEvSyn(ev)).norm4 = all.allEv(indxEvSyn(ev)).dist2centFSynch/maxdist2centFSynch;
        
    end
end