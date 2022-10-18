function [locEv, syn] = pieceOftrace2(locEv,roi, syn, InfoExp)
%% 20 before, 20 after
[~, firstPreframe] = max(InfoExp.descartFrame);
lastFrameBeforeTail = ((InfoExp.NumberAP)*((1/InfoExp.Frequency)/InfoExp.timeframe)) + firstPreframe;
lengthInterv = length(InfoExp.descartFrame); 
numFramBef = lengthInterv;
numFramAft = (lengthInterv*2) - 2;
numFramTot = numFramBef + numFramAft + 1;
evIndx = (1:1: length(locEv));



%% stuffLoc.allEv(Indexes_all_persyn(ev)).trace_dns = trace_dns(FrameN(ev)-numFramBef:(FrameN(ev) + numFramAft));
for synN = 1:length(roi)
    synID = roi(synN).synID;   
    traceDt = roi(synN).traceDt;
    mskSyn = ismember([locEv.synID], synID); 
    frameN = [locEv(mskSyn).frame];
    evIndxSyn = evIndx(mskSyn);
    mskSynch = [locEv(mskSyn).synch];
    mskAsynch2 = [locEv(mskSyn).asynch2nd];
    mskAsynchA = [locEv(mskSyn).asynchAlloc];
    
    avgTrace = cell(length(frameN),1);
    
    for ev = 1:length(frameN)
        if (frameN(ev)-numFramBef) > 0 && (frameN(ev)+numFramAft) <= lastFrameBeforeTail
            locEv(evIndxSyn(ev)).traceDt = (traceDt(frameN(ev)-numFramBef:frameN(ev)+numFramAft));
            avgTrace{ev} = (traceDt(frameN(ev)-numFramBef:frameN(ev)+numFramAft));
        else
            locEv(evIndxSyn(ev)).traceDt = NaN(numFramTot,1);
            avgTrace{ev} = NaN(numFramTot,1);
        end    
    end
  
    syn(synN).traceDtSynch =  nanmean(cell2mat((avgTrace(mskSynch))'),2);
    syn(synN).NumSynch = sum(mskSynch);
    %syn(synN).traceDtAsynch2 = nanmean(cell2mat(avgTrace(mskAsynch2)), 1);
    %syn(synN).traceDtAsynch2 = nanmean(cell2mat((avgTrace(mskAsynch2))'), 2);
    if sum(mskAsynch2) > 0
        %syn(synN).traceDtAsynch2 = nanmean(cell2mat(avgTrace(mskAsynch2)), 1);
        syn(synN).traceDtAsynch2 = nanmean(cell2mat((avgTrace(mskAsynch2))'), 2);
    else
        syn(synN).traceDtAsynch2 = NaN(numFramTot,1);
    end
    syn(synN).NumAsynch2 = sum(mskAsynch2);
    %syn(synN).traceDtAsynchA = nanmean(cell2mat(avgTrace(logical(mskAsynchA))), 1);
    syn(synN).traceDtAsynchA = nanmean(cell2mat((avgTrace(logical(mskAsynchA)))'), 2);
    syn(synN).NumAsynchAllo = sum(mskAsynchA);
    

end
%% replace [] for NaN

for ev1 = 1:length(locEv)
    if isempty(locEv(ev1).traceDt)
        locEv(ev1).traceDt = NaN;
    end        
end


end