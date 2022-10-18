function ALL = addDfNorm2ALL(ALL, synchFlag)
if synchFlag == 1
    evFr = 21;
else
    evFr = 22;
end

for i = 1 : length(ALL.allEv)
    
    trace = ALL.allEv(i).trace;
    trace = trace - min(trace) + 1;
    f0 = mean (trace(15:20));
    ALL.allEv(i).df = (trace(evFr) - f0)/ f0;
    ALL.allEv(i).trace_df = (trace - f0)/f0;
end