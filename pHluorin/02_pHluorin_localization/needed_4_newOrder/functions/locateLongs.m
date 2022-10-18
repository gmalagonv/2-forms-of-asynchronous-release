function  maskLongOnes = locateLongs (A)
%% Aug 2020, locate the events including more than one event. 
%% input: A, a cell with the amplitudes from fit
extended = [A{:}];
dimensions =  cellfun(@length, A);
lones = dimensions > 1;

maskLongOnes = false(length(extended),1); 
for i = 1:length(lones)
    if lones(i) == true
        if i == 1
            maskLongOnes(1:dimensions(i)) = true;
        else
            locPrev = sum(dimensions(1:(i-1)));
             maskLongOnes((locPrev + 1):(locPrev + dimensions(i)))= true;       
        end        
    end   
end
end