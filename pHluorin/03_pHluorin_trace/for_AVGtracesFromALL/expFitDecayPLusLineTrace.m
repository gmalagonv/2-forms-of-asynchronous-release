function [ratio, xunc, fitted, ySubsS, ySubsF] = expFitDecayPLusLineTrace(trace, prestim, lTrace, startPoint, filterFlag,subsFlag,varargin)
%t = (1:1:20)';
% if frameRelat == 0
%     lTrace = 20;
% else
%     lTrace = 21 - frameRelat;
% end
%% substract linear component
%subsFlag = 1; 
%startPoint = 5; % 0 used for ALL traces to highline the small difference in the fast component between synch and asynch.when 0, the statr point is sellected in each case




%% varargin, tau limits: varargin{1} = minTau, varargin{2} = maxTau 
%% numExp 1 or 2, single or double exponnential
t = (0:0.05:((lTrace -1)*0.05))';
y =  trace (prestim + 1:prestim + lTrace);

%% decide best startPoint
if startPoint == 0
    startPoints = [3,4,5,6,7,8, 9, 10];
    startPointsRes = zeros(size(startPoints));
    
    for stp = 1:length(startPoints)
        pnt = startPoints(stp);
        [slopeSTP, ~] = polyfit(t(pnt:end), (y(pnt:end)),1);
        ySTP= slopeSTP(1)*(t(pnt:end))+slopeSTP(2);
        startPointsRes(stp) = mean((y(pnt:end) - ySTP).^2);
        
    end
    [~, startPoint] =min(startPointsRes);
    startPoint = startPoints(startPoint);
end

%%
[slope, ~] = polyfit(t(startPoint:end), (y(startPoint:end)),1);
% 
% if abs(slope(1)) > 100
%    [slope, ~] = polyfit(t(5:end), (y(5:end)),1);
% %     varargin{1} = 0;
% %   filterFlag = 1;
% end

if filterFlag == 1
    yF = lowpass(y,0.5,20);
    [slope, ~] = polyfit(t(startPoint:end), yF(startPoint:end), 1);
end


yl= slope(1)*(t)+slope(2);
ySubsS = y - yl;
noise = 2*std(ySubsS(8:lTrace));


if subsFlag == 1
    slopeF = slope(1);
    intercF = slope(2);
    y = ySubsS;
    interc = 0;
    slope = 0;
   % disp('here')
    %% calculate min possible Tau
    yt = zeros(size(y));
    yt(1) = y(1);
    Ft = @(x,xdata)x(1)*exp(-x(2)*xdata);
    lbt = [yt(1)-noise -inf];
    ubt = [yt(1)+noise inf];
    if isempty(varargin) || varargin{1} == 0
        x0t = [yt(1) 0.1];
    else
        x0t = [yt(1) 1/varargin{1}];
    end
    options = optimoptions(@lsqcurvefit,'Display','none','MaxIterations',500);%,'UseParallel', true);
    [xunct, ~]= lsqcurvefit(Ft, x0t, t, yt, lbt , ubt , options);
    minPossTau = xunct(2);
    %% calculate max possible Tau
    yt = zeros(size(y));
    
    %% line decay
    [slm_ld, ~] = polyfit([1, startPoint], [y(1),0], 1);
    Ftt = @(xdata)xdata*slm_ld(1) + slm_ld(2);
    %%
    ld = Ftt(1:startPoint);
    yt(1:startPoint) = ld;
    [xunct, ~]= lsqcurvefit(Ft, x0t, t, yt, lbt , ubt , options);
    maxPossTau = xunct(2);
    
else
    interc = slope(2);
    slope = slope(1);
    
end

asymptote = mean(ySubsS(16:lTrace)) + interc;
%asymptote = mean(y(16:lTrace));

maxAsym = asymptote + noise;
minAsym = asymptote - noise;


%% test . By this there is no room for getting an asymptote during the fit of the exponential. 
if subsFlag == 1
    asymptote = 0;
    maxAsym = 0;
    minAsym = 0;
end
Ampli = y(1) - asymptote;%
%Ampli = max(y(1:4)) - asymptote;
maxAmpli = Ampli + noise;
minAmpli = Ampli - noise;






F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3) + x(4)*xdata;
%% VARIABLES
% x(1) = amplitude of tau1
% x(2) = 1 / tau1
% x(3) = asymptote
% x(4) = slope decay
%% minimim possible Tau: decays in 1 ponit from max to 0. 
%% in that case, t1/2 (time required for the decaying quantity to fall to one half of its initial value) 
%% t1/2  = 0.025 if sampling = 0.05. 
%% so, Tau = 0.025* (1/log(2)), so Tau = 0.0361; 



%%
if isempty(varargin) || varargin{1} == 0
    x0 = [Ampli 0.1 asymptote slope];
    if subsFlag == 0
        
        lb = [-inf -inf minAsym -inf];
        ub = [inf inf maxAsym inf];
        %ub = [inf 1/(0.025* (1/log(2))) maxAsym inf];
    else
%         lb = [-inf -inf minAsym slope];
%         ub = [inf inf maxAsym slope];
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       % lb = [maxAmpli 5 minAsym slope];
       % ub = [maxAmpli 500 maxAsym slope];
       %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       lb = [maxAmpli maxPossTau minAsym slope]; 
       ub = [maxAmpli minPossTau maxAsym slope];
        
    end



else
    Tau = 1/varargin{1};
    %slope = -1*varargin{2}; this worked with CT using 4 
    x0 = [Ampli Tau asymptote slope];
    
    if Tau > 0
        
        if subsFlag == 0
            lTau = Tau/4;
            uTau = Tau * 4;
        else
            lTau = maxPossTau;
            uTau = minPossTau;
        end
    else
        
        if subsFlag == 0
            lTau = 0;
            uTau = inf;
        else
            lTau = maxPossTau;
            uTau = minPossTau;
        end
    end
    
%% option A
%     if slope > 0
%         lSlope = -slope;
%         uSlope = slope;
%     else
%         lSlope = slope;
%         uSlope = -slope;
%     end
    
    %% option B thusi worked for CT with a factor of 4
%     sl(1) = ((yl(end) + noise/4) - (yl(6) - noise/4)) / (t(end) - t(6));
%     sl(2) = ((yl(end) - noise/4) - (yl(6) + noise/4)) / (t(end) - t(6));
%     
%     lSlope = min(sl);
%     uSlope = max(sl);
% %     
     %% option C
    lSlope = slope;
    uSlope = slope;
    
    lb = [minAmpli lTau minAsym lSlope];% interc];%*factSlopeUB];
    ub = [maxAmpli uTau maxAsym uSlope];% interc];%*(1/factSlopeUB)];

%         


end

%%

options = optimoptions(@lsqcurvefit,'Display','none','MaxIterations',500);%,'UseParallel', true);
[xunc, ~]= lsqcurvefit(F, x0, t, y, lb , ub , options);

if subsFlag == 0
    fitted =  F(xunc,t);
else
    fittedF = F(xunc,t);
   
    xuncF = xunc;
    xuncF(3) = xunc(3) + intercF;
    xuncF(4) = slopeF;
    fitted =  F(xuncF,t);
    y =  trace (prestim + 1:prestim + lTrace);
    ySubsF = y - fittedF;
    
     
    
end

%res = res/length(y);

%% normalize fitted and min sq


res = mean((y - fitted).^2);

minFitted = min(fitted);
fittedN = fitted - minFitted;
maxFitted = max(fittedN);
fittedN = fittedN ./maxFitted;

yN = y - minFitted;
yN = yN ./maxFitted;
resN = mean((yN - fittedN).^2);



% if fitted(1) > mean(fitted(2:3))
%     ratio = (fitted(1) - mean(fitted(2:3)))/ fitted(1);
% else
%     ratio = NaN;
% end

% if ratio > 0.8
%     disp('bigOne')
% end


%% reorder and ratio: first the fast tau a its ratio, then the slow one


%     xunc(1,1) = f_r;
xunc(1,2) = 1/xunc(1,2);
xunc(1,4) =-1*xunc(1,4);
xunc(1,5) = res;
xunc(1,6) = resN; 


if subsFlag == 1
    
    xunc(1,4) = -1*slopeF;
    xunc(1,5) = intercF;
    xunc(1,6) = fittedF(1);
%     if fittedF(1) < fittedF(2)
%         res = 1000000;
%     end
    xunc(1,7) = res;
    xunc(1,8) = resN; 

end


if subsFlag == 0  ||  xunc(1,4) < 0
    ratio = 0;
else 
    %fact = min([yl(1), fittedF(1)]);
    %ratio = (yl(1) -yl(end)) /(fittedF(1) - xunc(3)) ;%y(1);%
    %ratio = yl(1) /fittedF(1) ;%y(1);%
    ratio = yl(1) / fittedF(1) ;%y(1);%

end
if ratio < 0
    ratio = 0;
end
if ratio >= 10
    ratio = NaN;
end


% if ratio < 0 
%     disp('here')
% end
% if abs(ratio) > 20
%     ratio = NaN;
% 
% end
% clf(2), figure(2), plot(y),hold on, plot(fitted)%, plot((5:20),yF)
% disp('here')

% if fittedF(2) > fittedF(1) && startPoint <= 6
%     [ratio, xunc, fitted, ySubsS, ySubsF] = expFitDecayPLusLineTrace(trace, prestim, lTrace, startPoint +1, filterFlag,subsFlag,varargin{1});
% 
% end

end

%a = 1, figure(1), plot(t, [averages.bin.Traces(:,a)]), hold on, plot(t, F([averages.bin.fit(a,:)],t);
