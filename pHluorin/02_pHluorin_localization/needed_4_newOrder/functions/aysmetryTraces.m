function  [maxAsym, medAsym, varargout]= aysmetryTraces(xtrc, ytrc, tr1trc, tr2trc, includePeakFlag, fitFlag)
%% best corrType = 1;
corrType = 1;

if fitFlag == 1
    includePeakFlag = 0;
    corrType = 0; %% to keep all peaks at the same level
    %allInOne = mean([xtrc, ytrc, tr1trc, tr2trc],2);
end
    

%% trim borders;
borderL = 9; %11
trimmer = true(size(xtrc));
trimmer(1:borderL) = false;
trimmer(length(trimmer)-borderL+1:end) = false;

xtrc = xtrc(trimmer);
ytrc = ytrc(trimmer);
tr1trc = tr1trc(trimmer);
tr2trc = tr2trc(trimmer);


%halfLength = floor(length(xtrc)/2);
halfLoc = ceil(length(xtrc)/2);

if includePeakFlag == 1
    Lindex = (1:halfLoc);
    Rindex = ((halfLoc):length(xtrc));
else
    Lindex = (1:halfLoc -1);
    Rindex = ((halfLoc+1):length(xtrc));
    if fitFlag == 1
        peak = mean([xtrc, ytrc, tr1trc, tr2trc],2);
        peak = peak(halfLoc);
        
    end
    
end




%% parts ordered clock-wise
xL = xtrc(Lindex);
t2R = tr2trc(Rindex);
yR = ytrc(Rindex);
t1R = tr1trc(Rindex);
xR = xtrc(Rindex);
t2L = tr2trc(Lindex);
yL = ytrc(Lindex);
t1L = tr1trc(Lindex);
%ints = zeros(2,1);

%% total Axis = 4
axis = struct;

 %% Axis 1 = x
% axis(1).halfOneInt = trapz(t2R) + trapz(yR) + trapz(t1R);
% axis(1).halfTwoInt = trapz(t2L) + trapz(yL) + trapz(t1L);
if fitFlag == 0
    axis(1).tracehalfOne = (t2R + yR +t1R)/3;
    axis(1).tracehalfTwo = (t2L + yL +t1L)/3;

    axis(2).tracehalfOne = (yR + t1R + xR)/3;
    axis(2).tracehalfTwo = (yL + t1L + xL)/3;

    axis(3).tracehalfOne = (t1R + xR + flip(t2L))/3;
    axis(3).tracehalfTwo = (t1L + xL + flip(t2R))/3;

    axis(4).tracehalfOne = (xR + flip(t2L) + flip(yL))/3;
    axis(4).tracehalfTwo = (xL + flip(t2R) + flip(yR))/3;
else
    axis(1).tracehalfOne = [flip((t2R + yR +t1R)/3); peak ; ((t2R + yR +t1R)/3)];
    axis(1).tracehalfTwo = [(t2L + yL +t1L)/3 ; peak; flip((t2L + yL +t1L)/3)];

    axis(2).tracehalfOne = [flip((yR + t1R + xR)/3); peak; ((yR + t1R + xR)/3)];
    axis(2).tracehalfTwo = [((yL + t1L + xL)/3); peak; flip((yL + t1L + xL)/3)];

    axis(3).tracehalfOne = [flip((t1R + xR + flip(t2L))/3); peak; ((t1R + xR + flip(t2L))/3)];
    axis(3).tracehalfTwo = [((t1L + xL + flip(t2R))/3); peak; flip((t1L + xL + flip(t2R))/3)];

    axis(4).tracehalfOne = [flip((xR + flip(t2L) + flip(yL))/3); peak; ((xR + flip(t2L) + flip(yL))/3)];
    axis(4).tracehalfTwo = [((xL + flip(t2R) + flip(yR))/3); peak; flip((xL + flip(t2R) + flip(yR))/3)];
    
end

if fitFlag == 0
    for i = 1:4
    %     InthalfOne = trapz(axis(i).tracehalfOne);
    %     InthalfTwo = trapz(axis(i).tracehalfTwo);
        if corrType ~= 0
            if corrType == 1
                corr = min([min(axis(i).tracehalfOne), min(axis(i).tracehalfTwo)]);
                %corr = ones(size(axis(i).tracehalfOne)) .* factor;
                axis(i).tracehalfOne = axis(i).tracehalfOne - corr;
                axis(i).tracehalfTwo = axis(i).tracehalfTwo - corr;

            else
                ints(1,1) = trapz(axis(i).tracehalfOne);
                ints(1,2) = trapz(axis(i).tracehalfTwo);
                if sum(ints < 0) ~= 0
                    corr = ones(size(xR)).*((abs(min(ints)) + 1)/(length(xR)-1));
                    axis(i).tracehalfOne = axis(i).tracehalfOne + corr;
                    axis(i).tracehalfTwo = axis(i).tracehalfTwo + corr;
                end

            end
        end

        int0 = min(trapz(axis(i).tracehalfOne), trapz(axis(i).tracehalfTwo));
        int1 = max(trapz(axis(i).tracehalfOne), trapz(axis(i).tracehalfTwo));

        %axis(i).ratio = (int1 - int0)/ int0;

        if corrType == 1
            axis(i).ratio = (int1 - int0)/ int0;
        elseif corrType == 2
             axis(i).ratio = (int1 / int0);
        end
    end
else
    parfor  i = 1:4

        [~,fitParams1]= gaussFit3CoCenter(axis(i).tracehalfOne,  [0 halfLoc 0 0]);
        [~,fitParams2]= gaussFit3CoCenter(axis(i).tracehalfTwo,  [0 halfLoc 0 0]);

        axis(i).sigma0 = min(fitParams1(3),fitParams2(3));
        axis(i).sigma1 = max(fitParams1(3),fitParams2(3));
        axis(i).ratio = (axis(i).sigma1 / axis(i).sigma0);
        %axis(i).ratio = ((sigma1 - sigma0)/ sigma0);

    end
end  






[maxAsym, lc]= max([axis.ratio]);
medAsym = median([axis.ratio]);
if fitFlag == 1
    varargout{1} = axis(lc).sigma0;
    varargout{2} = axis(lc).sigma1;
end

% if sigma0 < 0.501 || sigma1 > 9.9
%     maxAsym = NaN;
%     medAsym = NaN;
%     sigma0 = NaN;
%     sigma1 = NaN;
% end

end


