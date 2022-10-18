function [locEv, syn] = distances2centers(locEv, syn, InfoExp)
%% oct 2020. distances and averages are calculated for every event regardless of their nature 
%%
synIDs = [syn.synID];%round(([syn.rowNum] .* [syn.columnNum])./([syn.sigma]));
evIndx = (1:1: length(locEv));
%% 

for synN = 1:length(synIDs)
    synID = synIDs(synN);   
    mskSyn = ismember([locEv.synID], synID); 
    evIndxSyn = (evIndx(mskSyn))';
    mskSynch = [locEv(mskSyn).synch];
    mskAsynch2 = [locEv(mskSyn).asynch2nd];
    mskAsynchA = [locEv(mskSyn).asynchAlloc];
    XY = [[locEv(mskSyn).x]' [locEv(mskSyn).y]'];    
    %% XY and indexes in a cell
    masked = cell(4,1);
    %% SYNCH
    masked{1,1} = XY(mskSynch, :);
    masked{1,2} = evIndxSyn(mskSynch, :);
    masked{1,3} = 'Synch';
    %% SYNCH + 2nd ASYNCH
    if sum(mskAsynch2) ~= 0
        masked{2,1} = XY((logical(mskAsynch2 + mskSynch)),:);
        masked{2,2} = evIndxSyn((logical(mskAsynch2 + mskSynch)),:);
        masked{2,3} = 'Synch2nd';
    else
        masked{2,1} = NaN;
        masked{2,2} = NaN;
        masked{2,3} = 'Synch2nd';
    end
        
    %% ASYNCH
    if sum(mskAsynch2 + mskAsynchA) ~= 0
        masked{3,1} = XY((logical(mskAsynch2 + mskAsynchA)),:);
        masked{3,2} = evIndxSyn((logical(mskAsynch2 + mskAsynchA)),:);
        masked{3,3} = 'Asynch';
    else
        masked{3,1} = NaN;
        masked{3,2} = NaN;
        masked{3,3} = 'Asynch';
    end
    %% ALL
    masked{4,1} = XY;
    masked{4,2} = evIndxSyn;
    masked{4,3} = 'All';
    
     if sum(mskAsynch2) ~= 0
        masked{5,1} = XY(mskAsynch2,:);
        masked{5,2} = evIndxSyn(mskAsynch2,:);
        masked{5,3} = 'Asynch2';
    else
        masked{5,1} = NaN;
        masked{5,2} = NaN;
        masked{5,3} = 'Asynch2';
    end
    
    
    
    for i = 1:4 % only 4 of the mask things, last one only asynch2nd
        if length(masked{i,1}) >= 2%~isnan(masked{i,1})
            %% PHYsiological center
            pairComparisonPhy(1,1:2) = mean(masked{i,1});
            syn(synN).(['centerPhyX'  masked{i,3}]) = pairComparisonPhy(1,1);
            syn(synN).(['centerPhyY'  masked{i,3}]) = pairComparisonPhy(1,2);
            %% 'Fisical center' and area        
            [pairComparisonF, area] = centroidHull(masked{i,1}(:,1),masked{i,1}(:,2));
            syn(synN).(['centerFX'  masked{i,3}]) = pairComparisonF(1,1);
            syn(synN).(['centerFY'  masked{i,3}]) = pairComparisonF(1,2);
            syn(synN).(['area'  masked{i,3}]) = area;
            %rootName = masked{i,3};
   
            %% mean distances to centerS calculated for EVERY event regardless its nature
            %meaNdist2centPhy = 0;
            %meaNdist2centF = 0;
            for j = 1:length(XY)
                %% PHY
                pairComparisonPhy(2,1) = XY(j,1);
                pairComparisonPhy(2,2) = XY(j,2);
                %locEv(masked{i,2}(j)).(['dist2centPhy'  masked{i,3}]) = pdist(pairComparisonPhy)*InfoExp.PixelSize;
                locEv(evIndxSyn(j)).(['dist2centPhy'  masked{i,3}]) = pdist(pairComparisonPhy)*InfoExp.PixelSize;
                %meaNdist2centPhy = meaNdist2centPhy + locEv(evIndxSyn(j)).(['dist2centPhy'  masked{i,3}]);
                %% F
                pairComparisonF(2,1) = XY(j,1);
                pairComparisonF(2,2) = XY(j,2);
                %locEv(masked{i,2}(j)).(['dist2centF'  masked{i,3}]) = pdist(pairComparisonF)*InfoExp.PixelSize;
                locEv(evIndxSyn(j)).(['dist2centF'  masked{i,3}]) = pdist(pairComparisonF)*InfoExp.PixelSize;
                %meaNdist2centF = meaNdist2centF + locEv(evIndxSyn(j)).(['dist2centF'  masked{i,3}]);        
            end
                %syn(synN).(['Mdist2centPhy'  masked{i,3}]) = meaNdist2centPhy / j;
                %syn(synN).(['Mdist2centF'  masked{i,3}]) = meaNdist2centF / j;

            
                
                %% distances to centerS calculated only for the partcular type
%             for j = 1:size(masked,1)
%                 meaNdist2centPhy = 0;
%                 meaNdist2centF = 0;
%                 XYm = masked{j,1};
%                 evIndxSynm = masked{j,2};
%                 for ev = 1: length(XYm)
%                     %% PHY
%                     pairComparisonPhy(2,1) = XYm(ev,1);
%                     pairComparisonPhy(2,2) = XYm(ev,2);
%                     locEv(evIndxSynm(ev)).([ 'dist2centPhy'  masked{i,3}]) = pdist(pairComparisonPhy)*InfoExp.PixelSize;
%                     meaNdist2centPhy = meaNdist2centPhy + locEv(evIndxSynm(ev)).(['dist2centPhy'  masked{i,3}]);
%                     %% F
%                     pairComparisonF(2,1) = XYm(ev,1);
%                     pairComparisonF(2,2) = XYm(ev,2);
%                     locEv(evIndxSynm(ev)).(['dist2centF'  masked{i,3}]) = pdist(pairComparisonF)*InfoExp.PixelSize;
%                     meaNdist2centF = meaNdist2centF + locEv(evIndxSynm(ev)).(['dist2centF'  masked{i,3}]);
%                 end
%                 syn(synN).([masked{j,3} 'Mdist2centPhy' rootName]) = meaNdist2centPhy / ev;
%                 syn(synN).([masked{j,3} 'Mdist2centF' rootName]) = meaNdist2centF / ev;
%                 
%             end



            
        else
            syn(synN).(['centerPhyX'  masked{i,3}]) = NaN;
            syn(synN).(['centerPhyY'  masked{i,3}]) = NaN;
            syn(synN).(['centerFX'  masked{i,3}]) = NaN;
            syn(synN).(['centerFY'  masked{i,3}]) = NaN;
            syn(synN).(['area'  masked{i,3}]) = NaN;
            for j = 1:length(XY)
                locEv(evIndxSyn(j)).(['dist2centPhy'  masked{i,3}]) = NaN;
                locEv(evIndxSyn(j)).(['dist2centF'  masked{i,3}]) = NaN;                        
            end
%             syn(synN).(['Mdist2centPhy'  masked{i,3}]) = NaN;
%             syn(synN).(['Mdist2centF'  masked{i,3}]) = NaN;           
        end
        
    end
    %% loop to get the avg distances per synapse, EXCLUSIVELY for Synch and Asynch2nd evs
    
    for k = 1:2
      
        syn(synN).(['m' masked{1,3} 'd2centFSynch']) = mean([locEv(masked{1,2}).dist2centFSynch]);
        if length(masked{5,1}) >= 2
            syn(synN).(['m' masked{5,3} 'd2centFSynch']) = mean([locEv(masked{5,2}).dist2centFSynch]);
        else
            syn(synN).(['m' masked{5,3} 'd2centFSynch']) = NaN;
        end

    end
end
end
    
    
