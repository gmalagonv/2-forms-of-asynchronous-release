function nuor_randSynch(Data, evType, dist2cRand, tethaRand, andorFlag, rmvAsynchAllo,iterations, varargin)%, frameRand, stimRand)
%% function to simulate evs
%% internal flag: how the ranges to get the random numbers for the dist2cent are calculated. 
% distRange == 1 or == 6 get the closer reults to the data
%distRange = 7;

%% by default distRange = 7. Otherwise set the number desired at the dist2cRand input. 
if length(dist2cRand) > 1    
    distRange = dist2cRand(2);
    dist2cRand = 1;
else
    distRange = 7;

end


count = 0;

% when synch events are simulated,  if == 1 a new center is calculated based on the hull envelope. Then the distances to that new center are calculated. 
% If the center is not replaced (replaceCenter = 0), all the locs are corected for both centers to match. In
% that case, the dist2cent of the asynch evs remain intact. 
replaceCenter = 0; % normally 0; 
removeFields = 1;%% remove innecesary fields to make smaller file

if iterations > 1
    removeFields = 1;
end



if andorFlag == 1
    pixSize = 86.666666666666670;
else
    pixSize = 65;
end

DataS = Data;
if rmvAsynchAllo == 1
    DataS.allEv = DataS.allEv(~[DataS.allEv.asynchAlloc]);
end

%% remove fields
if removeFields == 1
    Evfields = fieldnames(DataS.allEv);
    Evfields2keep = {'synID'; 'frame'; 'synch'; 'asynch2nd'; 'asynchAlloc'; 'x'; 'y'; 'dist2centFSynch'; 'tetha';'traceDt'};
    Evfields2rmv = Evfields(~(ismember(Evfields, Evfields2keep)));
    DataS.allEv = rmfield(DataS.allEv,Evfields2rmv);
    
    Synfields = fieldnames(DataS.perSyn);
    Synfields2keep = {'synID';'NumSynch'};
    Synfields2rmv = Synfields(~(ismember(Synfields, Synfields2keep)));
    DataS.perSyn = rmfield(DataS.perSyn,Synfields2rmv);
    if isfield(DataS, 'roi')
        DataS = rmfield(DataS, 'roi');
    end
end
%% several iterations (repeat several times the same subset of synapses)
if iterations > 1
    
    baseSyn = DataS.perSyn;
    baseEv = DataS.allEv;
 
    for it = 1:(iterations -1)
        syn2add = baseSyn;
        ev2add = baseEv;
        for z = 1:length(syn2add)
            newID =  str2double([num2str(syn2add(z).synID) num2str(it)]);
            syn2add(z).synID = newID;

        end
        syn2add = syn2add(~ismember([syn2add.synID],[DataS.perSyn.synID]));
        for w = 1:length(ev2add)
            ev2add(w).synID = str2double([num2str(ev2add(w).synID) num2str(it)]);
        end
        ev2add = ev2add(~ismember([ev2add.synID],[DataS.perSyn.synID]));
        
        
        DataS.perSyn = [DataS.perSyn syn2add];
        DataS.allEv = [DataS.allEv ev2add];
        
        
    end
end

indxAllEv =(1:length(DataS.allEv));

%% ALL events, used to get the distribution
if dist2cRand == 1
    if evType == 1
        distAllEvs2use = [Data.allEv([Data.allEv.synch]).dist2centFSynch];
    elseif evType == 2
        distAllEvs2use = [Data.allEv([Data.allEv.asynch2nd]).dist2centFSynch];
    else
        indxEvs2use = [Data.allEv(~[Data.allEv.synch]).dist2centFSynch];
    end
end




for syn = 1:length(DataS.perSyn)
    mskSyn = ismember([DataS.allEv.synID], DataS.perSyn(syn).synID);
    indxSynEv = indxAllEv(mskSyn);
    
    indxSynch = indxSynEv([DataS.allEv(indxSynEv).synch]);
    
    if evType == 1
        
        indxEvs2use = indxSynch;
        
    elseif evType == 2
        indxEvs2use = indxSynEv([DataS.allEv(indxSynEv).asynch2nd]);
    else
        indxEvs2use = indxSynEv(~[DataS.allEv(indxSynEv).synch]);
    end
    
    
    
    %% Center Location (needed to correct locations)
   if length(indxSynch) > 1 && length(indxEvs2use) >= 1
        xsynch = [DataS.allEv(indxSynch).x];
        ysynch = [DataS.allEv(indxSynch).y];
      

        indxHull = convhull(xsynch,ysynch);
        xHull = xsynch(indxHull);
        yHull = ysynch(indxHull);
        [centerF, ~] = centroidHull(xHull,yHull);


        %%%% polar coordinates

        if dist2cRand == 1
            if distRange < 7
                %% ranges of dist2c
                %% option1
                if distRange == 1
                    dist2cRange = [min([DataS.allEv(indxEvs2use).dist2centFSynch]) max([DataS.allEv(indxEvs2use).dist2centFSynch])];
                 %% option2
                elseif distRange == 2
                    dist2cRange = [0 max([DataS.allEv(indxEvs2use).dist2centFSynch])];
                %% option3
                elseif distRange == 3
                    meanDist2c= mean([DataS.allEv(indxEvs2use).dist2centFSynch]);
                    stdDist2c = std([DataS.allEv(indxEvs2use).dist2centFSynch]);
                    dist2cRange = [min([DataS.allEv(indxEvs2use).dist2centFSynch]) (meanDist2c + 2*stdDist2c)];
                %% option4
                elseif distRange == 4
                    meanDist2c= mean([DataS.allEv(indxEvs2use).dist2centFSynch]);
                    stdDist2c = std([DataS.allEv(indxEvs2use).dist2centFSynch]);
                    dist2cRange = [0 (meanDist2c + 2*stdDist2c)];        

                %% option5: use the dist2centFSynch from synch
                elseif distRange == 5
                    dist2cRange = [min([DataS.allEv(indxSynch).dist2centFSynch]) max([DataS.allEv(indxSynch).dist2centFSynch])];
                %% option6
                elseif distRange == 6
                 %dist2cRange = [min([DataS.allEv(indxEvs2use).dist2centFSynch])  1.5*max([DataS.allEv(indxSynch).dist2centFSynch])];
                    meanDist2c= mean([DataS.allEv(indxSynch).dist2centFSynch]);
                    stdDist2c = std([DataS.allEv(indxSynch).dist2centFSynch]);
                    %dist2cRange = [min([DataS.allEv(indxEvs2use).dist2centFSynch]) (meanDist2c + 2*stdDist2c)];
                    dist2cRange = [min([DataS.allEv(indxSynch).dist2centFSynch]) (meanDist2c + 2*stdDist2c)];   
                end
                dist2cent = ((dist2cRange(2) - dist2cRange(1)).*rand(length(indxEvs2use),1) + dist2cRange(1));
            else
                if distRange == 7
                    try
                        dist2cent = randFromdistribution(length(indxEvs2use), [DataS.allEv(indxEvs2use).dist2centFSynch]);
                    catch
                        dist2cent = randFromdistribution(length(indxEvs2use), distAllEvs2use);
                        count = count + 1;
                    end
                else
                    dist2cent = randFromdistribution(length(indxEvs2use), distAllEvs2use);                    
                end
                    
                
            end
            
        end
        if tethaRand == 1
            angles = 2.*rand(length(indxEvs2use),1)-1;
        end
        %% frame
        %% stim
        
       
        

        for ev = 1:length(indxEvs2use)
            evIndx = indxEvs2use(ev);
            if dist2cRand == 1
                dist = dist2cent(ev);
                if evType ~= 1
                    DataS.allEv(evIndx).dist2centFSynch = dist;
                end
                
            else
                dist = DataS.allEv(evIndx).dist2centFSynch;
                
            end

            if tethaRand == 1
                ang = angles(ev)* pi;
                DataS.allEv(evIndx).tetha = ang;
            else
                ang = DataS.allEv(evIndx).tetha;
            end
            [x,y] = pol2cart(ang,(dist/pixSize));

            DataS.allEv(evIndx).x = x + centerF(1);
            DataS.allEv(evIndx).y = y + centerF(2);          


        end
        
        %% recalculate centers, hull envelope, and dist2cents for ALL evs with the new coordinates if synch events are ramdomized
        if evType == 1
            xsynch = [DataS.allEv(indxSynch).x];
            ysynch = [DataS.allEv(indxSynch).y];
            indxHull = convhull(xsynch,ysynch);
            xHull = xsynch(indxHull);
            yHull = ysynch(indxHull);
            [centerF2, ~] = centroidHull(xHull,yHull);
            if replaceCenter == 1
                for ev = 1:length(indxSynEv)
                    evIndx = indxSynEv(ev);
                    DataS.allEv(evIndx).dist2centFSynch = pdist([centerF2;[DataS.allEv(evIndx).x DataS.allEv(evIndx).y]])*pixSize;
                end
            else
                centCorr = centerF2 - centerF;
                for ev = 1:length(indxSynch)
                    evIndx = indxSynch(ev);
                    DataS.allEv(evIndx).x = DataS.allEv(evIndx).x - centCorr(1);  
                    DataS.allEv(evIndx).y = DataS.allEv(evIndx).y - centCorr(2);
                end
                %% test
                
                xsynch = [DataS.allEv(indxSynch).x];
                ysynch = [DataS.allEv(indxSynch).y];
                indxHull = convhull(xsynch,ysynch);
                xHull = xsynch(indxHull);
                yHull = ysynch(indxHull);
                [centerF2, ~] = centroidHull(xHull,yHull);
                
                for ev = 1:length(indxSynEv)
                    evIndx = indxSynEv(ev);
                    DataS.allEv(evIndx).dist2centFSynch = pdist([centerF2;[DataS.allEv(evIndx).x DataS.allEv(evIndx).y]])*pixSize;
                end
            end
        end
   end
   
   
   
end
%% INFO simulation (GM, Aug 2021)
DataS.infoSimu.evType = evType;
DataS.infoSimu.dist2centFSynch = (dist2cRand == 1);
DataS.infoSimu.tetha = (tethaRand == 1);

if dist2cRand == 1 
    DataS.infoSimu.distRange = distRange;
end
if evType == 1
    DataS.infoSimu.replaceCenter = logical(replaceCenter);
end
%% do again normByRadius
DataS = normByRadius2(DataS, 1, pixSize);

%all = normByRadius2(all, extlPolygon, PixelSize)
%% do again Clustering
[DataS.perSyn, DataS.allEv] = allTypeOfClusters(DataS.perSyn, DataS.allEv,pixSize, 1, 1);
%%[syn, locEv] = allTypeOfClusters(syn, locEv,Data, reclustFlag, reclustFlagRand)

if distRange == 7 && count ~= 0
    disp(count/length(Data.perSyn))
end

if ~isempty(varargin)
    name = ['sim_' varargin{1}];
else
    name = 'sim';
end

assignin('base', name, DataS) 

end
