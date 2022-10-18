% to measure the duration of ON and OFF states

function [ONduration, OFFduration, indexesON]  = intervalsOnOff(traceIn)
indexesON = {};
ONduration = double.empty;
OFFduration = double.empty;
trsh = 0.5;%25;
boolON = 0;
boolOFF = 0;
intervalSampling = 1;%5e-2; %5e-05;
for i = 1:length(traceIn)
    if traceIn(i) >= trsh
        if boolON == 0
           ONduration(length(ONduration) + 1) = 1;
           indexesON{length(ONduration),1} = i; 
           
           boolON = 1;
           boolOFF = 0;
           
        else            
            ONduration(length(ONduration)) = ONduration(length(ONduration)) + 1;
            indexesON{length(ONduration),1} = [indexesON{length(ONduration), 1} i] ;
        end
    else
        if boolOFF == 0
            OFFduration(length(OFFduration)+1) = 1;
            boolOFF = 1;
            boolON = 0;            
        else
            OFFduration(length(OFFduration)) = OFFduration(length(OFFduration)) + 1;            
        end
        
    end
end
ONduration = ONduration .* intervalSampling;
OFFduration = OFFduration .* intervalSampling;

end