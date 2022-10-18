%SEM by GM, may 2019
function semval = sem_gm(wave, varargin)
%% sense = 1 : per line; sense = 2 : per column
% remove NaNs
[a, b] = size(wave);
if ~isempty(varargin)
    sense = varargin{1};
    

    if sense == 1
        numCycles = a;
    else
        numCycles = b;
        wave = wave';
    end
end

if a > 1 && b > 1
    semval = zeros(numCycles,1);
    for i = 1:numCycles
        waveTemp = wave(i,:);
        waveTemp = waveTemp(~isnan(waveTemp));
        semval(i,1) = std(waveTemp)/sqrt(length(waveTemp));
        
    end
else
    wave = wave(~isnan(wave));
    semval = std(wave)/sqrt(length(wave));
end
if ~isempty(varargin) &&  sense == 2
    semval = semval';

end
end

