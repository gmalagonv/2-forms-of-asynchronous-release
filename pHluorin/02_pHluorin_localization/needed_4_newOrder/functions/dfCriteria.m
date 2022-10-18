function locEv = dfCriteria (locEv, roi, InfoExp,flagDepurate, flagCrit)
%% remove evs before the first stim
%[~, firstPreframe] = max(InfoExp.descartFrame);

firstPreframe = 4;
framesPerStim = (1/InfoExp.Frequency)/InfoExp.timeframe;
locEv = locEv([locEv.frame] > firstPreframe);

synID = unique([locEv.synID]', 'stable');
dfALL = zeros(size(locEv))';
roiIndx = (1:1:length(roi))';


for syn = 1: length(synID)
    mskSyn = ismember([locEv.synID], synID(syn))';
    framesSyn = [locEv(mskSyn).frame]';
    stim = floor((framesSyn - (firstPreframe +1))/framesPerStim) + 1;
    preStim = (stim - 1)*framesPerStim + firstPreframe;
    preframesSynStart = preStim -(firstPreframe-1);    
    preframesSynEnd = preStim;
    
    f0Syn = zeros(length(preframesSynStart),1);
    
    roimsk = ismember([roi.synID],synID(syn));
    roiIndxSyn = roiIndx(roimsk);
    
    trace = roi(roiIndxSyn).traceDt;
    %trace = trace - min(trace) + 1; %% only positive values > 0

    for j = 1:length(preframesSynStart)
        f0Syn(j,1) = mean(trace(preframesSynStart(j):preframesSynEnd(j)));    
    end

    fSyn = trace(framesSyn);
    df = (fSyn-f0Syn)./f0Syn;
    dfALL(mskSyn) = df;
    
    

end



%% remove the negative values
% negDFloc = [locEv.df] <= 0;
% negDFloc =  std([[locEv(negDFloc).df],[locEv(negDFloc).df]*-1]);
% accDFcrit = [locEv.df] > 2* negDFloc;
%%
negDFloc = dfALL <= 0;
negDFloc = std([dfALL(negDFloc),(dfALL(negDFloc)*-1)]);
accDFcrit = dfALL > 2* negDFloc;

%% add df

for i = 1:length(dfALL)
    locEv(i).df = dfALL(i);
     locEv(i).dfCrit = accDFcrit(i);
end
%%


if flagDepurate == 1 
    if flagCrit == 1
        mskPositdf = dfALL > 0;
        
    else
        mskPositdf = accDFcrit;
    end
    locEv = locEv(mskPositdf);
end




% a = [locEvC.frame]';
% b = (a-(floor(a/20))*20);
% b(b == 0) = 20;


end 