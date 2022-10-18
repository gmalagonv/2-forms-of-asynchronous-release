%function [duration,mskSyncAsync, mskSameLoc,durationCorr,dist, rejec]  = clasifyEvents(locEv,InfoExp)
function locEv  = clasifyEvents(locEv,InfoExp)
framesALL = [locEv.frame]';

%% check continuity & duration of every event

duration = ones(size(framesALL))';

framesNext_i = framesALL;
framesNext_i(1) = [];
framesNext_i(end + 1) = NaN;
framesNext_i = framesNext_i - framesALL;

framesPrev_i= framesALL;
framesPrev_i= [NaN ; framesPrev_i];
framesPrev_i(end) = [];
framesPrev_i = framesALL - framesPrev_i;


mskContinuity = logical((framesPrev_i == 1) + (framesNext_i == 1));
mskContinuityB = (logical(framesPrev_i == 1)  + logical(framesNext_i==1)) == true;
[ONduration, ~, ONindexes]  = intervalsOnOff(mskContinuity);


for i = 1:length(ONduration)
    numInterv = sum(mskContinuityB(min(ONindexes{i}):max(ONindexes{i})))/2;
    if numInterv == 1
        for j = 1:ONduration(i)
            if j == 1
                duration(ONindexes{i}(1)) = ONduration(i);
            else
                duration(ONindexes{i}(j)) = 0;
                
            end
            
        end        
    else
       firstLast = ONindexes{i}(mskContinuityB(min(ONindexes{i}):max(ONindexes{i})));
       first = firstLast(1:2:end);
       dur = (firstLast(2:2:end)-first);        
        for k = 1: numInterv
            duration(first(k)) = dur(k) + 1;            
            for l = 1:dur(k)
                duration(first(k) + l) = 0;
            end
        end
        
    end

end
duration = duration';

%% check locations for the frames inside each continuous
% mskLongDur = duration > 1;
% indxLongDur = (indx(mskLongDur))';
 mskSameLoc = false(size(duration));
%  msk2remove = true(size(duration));
% mskSameLoc(indxLongDur) = true;
% 
% for m = 1:length(indxLongDur)
%     for n = 1:(duration(indxLongDur) -1)
%     end
%      
% end
%  
dist = NaN(size(duration));
rejec = dist;
%limitDist = 50/InfoExp.PixelSize;
durationCorr = duration;

for m = 1: length(duration)
    if duration(m) > 1
        mskSameLoc(m) = true;
        rejected = duration(m) - 1;
        for fr = 1:(duration(m) - 1)            
            Indx2Comp = m + fr;
            IndxStand = Indx2Comp - 1;
            xy_pstd = sqrt(((locEv(IndxStand).x_pstd).^2)+((locEv(IndxStand).y_pstd).^2));            
            distance = pdist_several(locEv(IndxStand).x, locEv(IndxStand).y, locEv(Indx2Comp).x, locEv(Indx2Comp).y);
            dist(Indx2Comp) =  distance(3);
            
            if distance(3) <= xy_pstd %limitDist %
                mskSameLoc(Indx2Comp) = true;
                rejected = rejected - 1; 
            else
%                 for p = Indx2Comp:Indx2Comp + rejected - 1
%                     msk2remove(p) = false;                    
%                 end
                break
            end           
        end
        rejec(m) = rejected;
        durationCorr(m) = durationCorr(m) - rejected;
        
    end
end
msk2remove = logical((duration == 1) + mskSameLoc);

%% remove events
locEv = locEv(msk2remove);
durationCorr = durationCorr(msk2remove);
framesALL = framesALL(msk2remove);

%% SYNCHRONOUS EVENTS


indx = (1:1:length(framesALL));
[~, firstPreframe] = max(InfoExp.descartFrame);
NumframesperStim = (1/InfoExp.Frequency)/InfoExp.timeframe;
%%

framesSynch = (firstPreframe:NumframesperStim:(NumframesperStim*(InfoExp.NumberAP-1) + firstPreframe))' + 1;

maskSync = ismember(framesALL, framesSynch);
indxSync = indx(maskSync)';

%FrmsBefSync = framesALL((indxSync-1));
FrmsBefSync = (indxSync-1);
FrmsBefSync(FrmsBefSync == 0) = 1;
FrmsBefSync = framesALL(FrmsBefSync);

ExpectFrmBefSync = (framesALL(maskSync))-1;
maskBefSync = (FrmsBefSync-ExpectFrmBefSync) ~= 0;

indxAccSync = indxSync(maskBefSync);
finalMskSync = false(size(indx))';
finalMskSync(indxAccSync) =  true;


%% ASYNCHRONOUS EVENTS 2nd frame
framesAsynch = (firstPreframe:NumframesperStim:(NumframesperStim*(InfoExp.NumberAP-1) + firstPreframe))' + 2;

maskAsync = ismember(framesALL, framesAsynch);
indxAsync = indx(maskAsync)';

%FrmsBefAsync = framesALL(indxAsync-1);
FrmsBefAsync = indxAsync-1;
FrmsBefAsync(FrmsBefAsync == 0) = 1;
FrmsBefAsync = framesALL(FrmsBefAsync);

ExpectFrmBefAsync = (framesALL(maskAsync))-1;
maskBefAsync = (FrmsBefAsync-ExpectFrmBefAsync) ~= 0;
indxAccAsync = indxAsync(maskBefAsync);
finalMskAsync = false(size(indx))';
finalMskAsync(indxAccAsync) =  true;

mskSyncAsync = logical(finalMskSync + finalMskAsync);
mskAlllocAsync = logical((durationCorr ~= 0) .* (~mskSyncAsync));

%

for y = 1: length(locEv)
    locEv(y).synch = finalMskSync(y);
    locEv(y).asynch2nd = finalMskAsync(y);
    locEv(y).asynchAlloc = mskAlllocAsync(y);
    locEv(y).permanence = durationCorr(y);

end
 %% remove the subsequent events (permanence = 0) they come after a ~= 1 permanence frames
mskPermanence = (durationCorr ~= 0);
 locEv = locEv(mskPermanence);
 

end