function locEv = bestFitLongs(locEv)
%% when there are multiple localizations and fits in a single frame, long flag is true. 
%% In these cases, only one fit & localization is selcted using this function.
%%
totNumEv = length(locEv);
maskTot = true(totNumEv,1);
indx = (1:1:totNumEv);
originalLong = [locEv.long]'; 


%% generate frameID by concatenating synID & frame (GM, 10.09.2020)
%frameID = ([locEv.synID] .* [locEv.frame] ./ ([locEv.frameRelat]));
frameID = frameIDgenerator(locEv);
frameLongID = unique(frameID([locEv.long]), 'stable');

%disp (num2str(length(frameLongID)));
%kepti = cell(length(frameLongID),1);

for i = 1:length(frameLongID)
    localindx = indx(frameLongID(i) == frameID );
    if sum(diff(localindx)) == (length(localindx) -1) % only consequtive events
        %localindx = indx(ismember(frameID, frameLongID));
        %kepti{i} = localindx;
        localmask = false(size(localindx));
        localfactor = ([locEv(localindx).x_pstd] + [locEv(localindx).A_pstd]) + (1./[locEv(localindx).A]);
        [~, selectedfactor] = min(localfactor);
        %kepti(i) = length(selectedfactor);
        localmask(selectedfactor) = true;    
        maskTot(localindx(~localmask)) = false;    
        locEv(localindx(localmask)).long = false;
    else
        disp(['i = ' num2str(i) ',' ' something weird happened, redefine ID frames'])
        break
    end
    
end

locEv = locEv(maskTot);
if sum([locEv.long]) == 0 && sum(originalLong(maskTot)) == length(frameLongID)
    disp ('selected the best fit between multiple locs.');
    locEv = rmfield(locEv, 'long');
else
    %disp('something weird happened, check ID frames')
end
end