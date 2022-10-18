function [syn, locEv] = allTypeOfClusters(syn, locEv,Data, reclustFlag, reclustFlagRand)
%% reclustFlagRand if == 1, as many synch as asynch there are, are selected for recluster. They are selected randomly. 
if ~isstruct(Data)
    pixSize = Data;
    Data = struct;
    Data.InfoExp.PixelSize = pixSize;
    Data.InfoExp.timeframe = 0.05;
    
end
a = 0;   



%%
synIDs = [syn.synID];%round(([syn.rowNum] .* [syn.columnNum])./([syn.sigma]));
evIndx = (1:1: length(locEv));
%% 
%fields2rm = {'loc', 'maxRad', 'frames'};


for synN = 1:length(synIDs)
    synID = synIDs(synN);   
    mskSyn = ismember([locEv.synID], synID); 
    evIndxSyn = evIndx(mskSyn);
    mskSynch = [locEv(mskSyn).synch];
    %evIndxSynSynch = evIndxSyn(mskSynch);
    mskAsynch2 = [locEv(mskSyn).asynch2nd];
    mskAsynchA = [locEv(mskSyn).asynchAlloc];    
    XY = [[locEv(mskSyn).x]' [locEv(mskSyn).y]'];    
    frames = [locEv(mskSyn).frame]';
    
%% synch clusters
    if sum(mskSynch) > 2
        framesSynch = frames(mskSynch);
        XYsynch = XY(mskSynch, :);
        clustSynch = ClusterVesicFusion( Data, 0, 50, 0, 0, XYsynch);
        %% add the frame (time) when synch cluster was used. GM, may 2021, BOG
        clustSynch.cluster(1).frames = [];
        for evXY = 1:length(XYsynch)
            cl = 1;
            while ~ismember(XYsynch(evXY),  [clustSynch.cluster(cl).loc])
                cl = cl + 1;
            end
            clustSynch.cluster(cl).frames = [clustSynch.cluster(cl).frames framesSynch(evXY)];
        end
        
        
        syn(synN).clustSynch = clustSynch.cluster;
        syn(synN).NumclustSynch = length(syn(synN).clustSynch);
        
        centerSynch= [syn(synN).clustSynch.centroid];
        centerSynch = [centerSynch(1:2:length(centerSynch))',centerSynch(2:2:length(centerSynch))'];
        
       
        
        %centerSynch = cell2mat(squeeze(struct2cell(rmfield((syn(synN).clustSynch),fields2rm))));
    else
        syn(synN).clustSynch = NaN;
        syn(synN).NumclustSynch = NaN;
        centerSynch = NaN;
        
    end
%% synch & Asynch2nd clusters    
    if sum(mskAsynch2) > 0
        %XYasynch2nd = XY(mskAsynch2, :);
        clustSynch2nd = ClusterVesicFusion( Data, 0, 50, 0, 0, XY((logical(mskAsynch2 + mskSynch)),:));
        syn(synN).clustSynch2nd = clustSynch2nd.cluster;
        syn(synN).NumclustSynch2nd = length(syn(synN).clustSynch2nd);
        
        centerSynch2nd= [syn(synN).clustSynch2nd.centroid];
        centerSynch2nd = [centerSynch2nd(1:2:length(centerSynch2nd))',centerSynch2nd(2:2:length(centerSynch2nd))'];
        
        %centerSynch2nd = cell2mat(squeeze(struct2cell(rmfield((syn(synN).clustSynch2nd),fields2rm))));
%% REcluster removing in each Syn as many events as 2CC events are
%% to then calculate the distance from the removed events to the centers of the new clusters.   
        if reclustFlag == 1 && reclustFlagRand == 0
            centerReclusts = cell(length(XYsynch), 1);
            
            for i = 1: length(XYsynch)
                XYrem = XYsynch;
                XYrem(i, :) = [];
                REClust = ClusterVesicFusion( Data, 0, 50, 0, 0, XYrem);
                
                REClust= [REClust.cluster.centroid];
                centerReclusts{i,1}  = [REClust(1:2:length(REClust))',REClust(2:2:length(REClust))'];
            end
              %% option2
%             if reclustFlagRand == 1
%                 if sum(mskSynch) > sum(mskAsynch2)
%                     num2reclust =  randperm(sum(mskSynch));
%                     num2reclust = sort(num2reclust(1:sum(mskAsynch2)));
%                 else
%                     num2reclust = 1:sum(mskSynch);
%                 end
%                 
% 
%             end
 
       %% option1     
        elseif reclustFlag == 1 && reclustFlagRand == 1 

             
             if sum(mskSynch) > sum(mskAsynch2)
                 num2reclust =  randperm(sum(mskSynch));
                 num2reclust = sort(num2reclust(1:sum(mskAsynch2)));
             else
                 num2reclust = 1:sum(mskSynch);
             end
             centerReclusts = cell(length(num2reclust), 1);
             for i = 1:length(num2reclust)
                 XYrem = XYsynch;
                 XYrem(num2reclust(i), :) = [];
                 REClust = ClusterVesicFusion( Data, 0, 50, 0, 0, XYrem);
                 REClust= [REClust.cluster.centroid];
                 centerReclusts{i,1}  = [REClust(1:2:length(REClust))',REClust(2:2:length(REClust))'];

             end
            
        end        
        
    else 
        syn(synN).clustSynch2nd = NaN;
        syn(synN).NumclustSynch2nd = NaN;
        centerSynch2nd = [NaN NaN];
        %centerReclusts = {NaN NaN};
    end

%% Asynch clusters
    if sum(mskAsynchA) > 0
        clustAsynch = ClusterVesicFusion( Data, 0, 50, 0, 0, XY((logical(mskAsynch2 + mskAsynchA)),:));
        syn(synN).clustAsynch = clustAsynch.cluster;
        syn(synN).NumclustAsynch = length(syn(synN).clustAsynch);
        
        centerAsynch= [syn(synN).clustAsynch.centroid];
        centerAsynch = [centerAsynch(1:2:length(centerAsynch))',centerAsynch(2:2:length(centerAsynch))'];
        
        %centerAsynch = cell2mat(squeeze(struct2cell(rmfield((syn(synN).clustAsynch),fields2rm))));
    else
        syn(synN).clustAsynch = NaN;
        syn(synN).NumclustAsynch = NaN;
    end
%% ALL clusters
    if length(XY) > 2
        clustAll  = ClusterVesicFusion( Data, 0, 50, 0, 0, XY);
        syn(synN).clustAll = clustAll.cluster;
        syn(synN).NumclustAll = length(syn(synN).clustAll);
        
        centerAll = [syn(synN).clustAll.centroid];
        centerAll = [centerAll(1:2:length(centerAll))',centerAll(2:2:length(centerAll))'];
        %centerAll = cell2mat(squeeze(struct2cell(rmfield((syn(synN).clustAll),fields2rm))));
    else
        syn(synN).clustAll = NaN;
        syn(synN).NumclustAll = NaN;
        
    end
    
    %% ADD dist dist2ClustCent to each event   
    reclustIndx = 1;
    for ev = 1:length(evIndxSyn)
        %% synch clusters  & recluster
        %% distance and time since use of the closest synch cluster that is different to the one where the ev occurred (if not, for all synch == 0)
        if length(XY) > 2
            

            %% time and distance since last synch cluster used before & after ev
            
            diffFrames = [([syn(synN).clustSynch.frames])' (zeros(length([syn(synN).clustSynch.frames]),1))];
            
            for mm = 1:length(syn(synN).clustSynch)
                diffFrames((ismember (diffFrames(:,1), syn(synN).clustSynch(mm).frames)),2) = pdist ([centerSynch(mm,:); XY(ev,:)]) * Data.InfoExp.PixelSize;%centerSynch(mm,1);
            end
            %% remove cluster used by itself (only true for synch events)
            
            %% <-- test
           diffFrames(:,1) = diffFrames(:,1) - frames(ev);

%             if ismember(XY(ev),XYasynch2nd)
%                 diffFrames(:,1) = diffFrames(:,1) - (frames(ev) - 1);
%             else
%                 diffFrames(:,1) = diffFrames(:,1) - frames(ev);
%             end
            
            diffFrames = diffFrames((diffFrames(:,1) ~= 0), :);
            
            %%  OVERALL SPATIALLY CLOSEST CLUSTER USED
            if ismember(XY(ev), XYsynch)
               diffFramesCl = diffFrames((diffFrames(:,2) > 40), :);
            else
                diffFramesCl = diffFrames;
            end
            [DistClosestClust, ~] = min(diffFramesCl(:,2));
              
            locEv(evIndxSyn(ev)).dist2closestClust = DistClosestClust;
              
            mskDF = ismember(diffFramesCl(:,2), DistClosestClust);
            locEv(evIndxSyn(ev)).time2closestClust = min(abs(diffFramesCl(mskDF,1)));
            %% add fields to keep order
            locEv(evIndxSyn(ev)).dist2closestClustBef = [];
            locEv(evIndxSyn(ev)).time2closestClustBef = [];
            locEv(evIndxSyn(ev)).dist2closestClustAft = [];
            locEv(evIndxSyn(ev)).time2closestClustAft = [];
            
            %% OVERALL TEMPORALLY CLOSEST CLUSTER USED
            [valMin, locMin] = min(abs(diffFrames(:,1)));
            locEv(evIndxSyn(ev)).time2tempClosestClust = valMin * Data.InfoExp.timeframe;
            locEv(evIndxSyn(ev)).dist2tempClosestClust = diffFrames(locMin,2);
            %% add fields to keep order
            locEv(evIndxSyn(ev)).time2ClustBef = [];
            locEv(evIndxSyn(ev)).dist2ClustBef = [];
            locEv(evIndxSyn(ev)).time2ClustAft = [];
            locEv(evIndxSyn(ev)).dist2ClustAft = [];


            %% EVENTS BEFORE 
            mskBefore = diffFrames(:,1) < 0;
            if sum(mskBefore) ~= 0
                evsBefore = diffFrames(mskBefore,:);
                %% closest cluster (in distance) used before
                
                if ismember(XY(ev), XYsynch)
                    diffFramesCl = evsBefore((evsBefore(:,2) > 40), :);
                else
                    diffFramesCl = evsBefore;
                end
        

                [DistClosestClust, ~] = min(diffFramesCl(:,2));

                locEv(evIndxSyn(ev)).dist2closestClustBef = DistClosestClust;

                mskDF = ismember(diffFramesCl(:,2), DistClosestClust);
                locEv(evIndxSyn(ev)).time2closestClustBef = min(abs(diffFramesCl(mskDF,1)));
                
         
                %% time 2 cluster used immediatly before
                [valMinBefore, locMiBefore] = min(abs(evsBefore(:,1)));
                locEv(evIndxSyn(ev)).time2ClustBef = valMinBefore * Data.InfoExp.timeframe;
                locEv(evIndxSyn(ev)).dist2ClustBef = evsBefore(locMiBefore,2);
      
            end
            
            
            %% EVENTS AFTER
            mskAfter = diffFrames(:,1) > 0;
            if sum(mskAfter) ~= 0
                evsAfter = diffFrames(mskAfter,:);
                %% closest cluster used After
                if ismember(XY(ev), XYsynch)
                    diffFramesCl = evsAfter((evsAfter(:,2) > 40), :);
                else
                    diffFramesCl = evsAfter;
                end
                [DistClosestClust, ~] = min(diffFramesCl(:,2));

                locEv(evIndxSyn(ev)).dist2closestClustAft = DistClosestClust;

                mskDF = ismember(diffFramesCl(:,2), DistClosestClust);
                locEv(evIndxSyn(ev)).time2closestClustAft = min(abs(diffFramesCl(mskDF,1)));
                %% time 2 cluster used immediatly after
                
                [valMinAfter, locMiAfter] = min(abs(evsAfter(:,1)));
                locEv(evIndxSyn(ev)).time2ClustAft = valMinAfter * Data.InfoExp.timeframe;
                locEv(evIndxSyn(ev)).dist2ClustAft = evsAfter(locMiAfter,2);

            end
            
    
        else
            locEv(evIndxSyn(ev)).dist2ClustSynch = NaN;
        end
%% RECLUSTER
        if reclustFlag == 1 && sum(mskAsynch2) > 0
            if mskSynch(ev) == true %|| mskAsynch2(ev) == true
                if reclustFlagRand == 0
                    cents = centerReclusts{reclustIndx,1};
                    if ismember(XY(ev,1), cents(:,1))
                        disp('error with reclustering')
                    end
                    dist = pdist_several(XY(ev,1),XY(ev,2),cents(:,1),cents(:,2));
                    locEv(evIndxSyn(ev)).dist2Reclust = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
                    reclustIndx = reclustIndx + 1;
                    
                elseif reclustFlagRand == 1 && reclustIndx <= length(num2reclust) && XY(ev) == XYsynch(num2reclust(reclustIndx))
                    cents = centerReclusts{reclustIndx,1};
                    if ismember(XY(ev,1), cents(:,1))
                        disp('error with reclustering')
                    end
                    dist = pdist_several(XY(ev,1),XY(ev,2),cents(:,1),cents(:,2));
                    locEv(evIndxSyn(ev)).dist2Reclust = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
                    reclustIndx = reclustIndx + 1;
                    
                    
                
                else
                    locEv(evIndxSyn(ev)).dist2Reclust = NaN;

                end
            elseif mskAsynch2(ev) == true
                %% for asynch events, select randomly one of the cluster configurations.
                [numReclusts, ~] = size(centerReclusts);
                cent2use = randi(numReclusts);
                cents = centerReclusts{cent2use,1};
                dist = pdist_several(XY(ev,1),XY(ev,2),cents(:,1),cents(:,2));
                locEv(evIndxSyn(ev)).dist2Reclust = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
            else
                locEv(evIndxSyn(ev)).dist2Reclust = NaN;
            end
        else
            locEv(evIndxSyn(ev)).dist2Reclust = NaN;

        end
%% END RECLUSTER

        %%  Asynch2nd clusters %% check cases when cluaster is 'mixed'. what proportion per synapse?
        if mskAsynch2(ev) > 0 || mskSynch(ev) > 0                         
            dist = pdist_several(XY(ev,1),XY(ev,2),centerSynch2nd(:,1),centerSynch2nd(:,2));
            %[~ , mnLoc]= min(dist(:,3));
            dist = dist(~(dist(:,3) == min(dist(:,3))), :);
            %dist = dist((dist(:,3) ~= 0),:);
            
            locEv(evIndxSyn(ev)).dist2ClustSynch2nd = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
        else
            locEv(evIndxSyn(ev)).dist2ClustSynch2nd = NaN;
            
        end
        %% Asynch clusters
         if sum(mskAsynchA) > 0
            dist = pdist_several(XY(ev,1),XY(ev,2),centerAsynch(:,1),centerAsynch(:,2));
            locEv(evIndxSyn(ev)).dist2ClustAsynch = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
         end
        %% ALL clusters
        dist = pdist_several(XY(ev,1),XY(ev,2),centerAll(:,1),centerAll(:,2));
        locEv(evIndxSyn(ev)).dist2ClustAll = (min(dist(:,3)))*(Data.InfoExp.PixelSize);
    
    
    end 

    %% add to syn a field with the mean Reclust dist per synapse, both synch and asynch
    if reclustFlag == 1
        syn(synN).meanDistReclSynch = mean([locEv(evIndxSyn(mskSynch)).dist2Reclust], 'omitnan');
        syn(synN).meanDistReclAsynch = mean([locEv(evIndxSyn(mskAsynch2)).dist2Reclust], 'omitnan');
    end
 
end
end
