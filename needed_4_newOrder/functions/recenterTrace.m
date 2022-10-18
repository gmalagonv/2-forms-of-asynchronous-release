function trace = recenterTrace(trace, locCenter, idealLocCenter)
sizeBorders = 3;
lengTrace = length(trace);
mskBorders = false(size(trace));
mskBorders(1:sizeBorders) = true;
mskBorders((lengTrace-sizeBorders):end) = true;
replacer = mean(trace(mskBorders)) * ones(size(trace));

diff = idealLocCenter - locCenter;

if diff < 0 %% it needs to be moved to the right
    diff = abs(diff);
    replacer(1:(lengTrace - diff)) = trace((1 + diff):end);
    
else
    diff = abs(diff);
    replacer((1 + diff):end) = trace(1:(lengTrace - diff));
 
    
end
trace = replacer;



end
