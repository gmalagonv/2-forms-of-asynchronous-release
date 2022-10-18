function Data = transTraceInt(DataName,rmvBckgFlag,foldFlag, noiseFlag)

if ischar(DataName)
    Data = evalin('base', DataName);
else
    Data = DataName;
end

if length([Data.perSyn.synID])== length(unique([Data.perSyn.synID])) && length([Data.perSyn.synID]) == length(unique([Data.allEv.synID]))
else
    Data = corrDuplisynIDS(Data);
end



for synNum = 1:length(Data.roi)
    synID = Data.roi(synNum).synID;
    mskSyn = [Data.allEv.synID] == synID;
    movie = Data.roi(synNum).movie;
    if synNum == 1
        locEv = 0;

    end
    
    for ev = 1:sum(mskSyn)
        frameN = Data.allEv(locEv + ev).frame;
        x = round(Data.allEv(locEv + ev).x);
        y = round(Data.allEv(locEv + ev).y);
        
        
        if noiseFlag == 1
            x = x + 10;
            y = y + 10;
            if x > 41
                x = x - 20;
            end
            
            if y > 41
                y = y - 20;
            end
        end
        if noiseFlag == 0
            Data.allEv(locEv + ev).traceInt  = traceInt(movie, x, y,frameN,rmvBckgFlag,foldFlag, 20, 38);
            locPeak = 21;
            Data.allEv(locEv + ev).traceIntPeak = Data.allEv(locEv + ev).traceInt(locPeak) - mean(Data.allEv(locEv + ev).traceInt((locPeak-6):(locPeak-1)));
            Data.allEv(locEv + ev).traceIntrat1 = Data.allEv(locEv + ev).traceInt(locPeak)/Data.allEv(locEv + ev).traceInt(locPeak + 1);
            Data.allEv(locEv + ev).traceIntrat2 = Data.allEv(locEv + ev).traceInt(locPeak + 1)/Data.allEv(locEv + ev).traceInt(locPeak + 2);
            
            
            
        else
            Data.allEv(locEv + ev).traceIntN  = traceInt(movie, x, y,frameN,rmvBckgFlag,foldFlag, 20, 38);
            
        end
        
    end
    locEv = locEv + ev;
end
    