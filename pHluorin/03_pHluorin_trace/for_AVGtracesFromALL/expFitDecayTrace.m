function [ratio, xunc, fitted] = expFitDecayTrace(trace, prestim, numExp,varargin)
%t = (1:1:20)';
%% varargin, tau limits: varargin{1} = minTau, varargin{2} = maxTau 
%% numExp 1 or 2, single or double exponnential
t = (0:0.05:(19*0.05))';
y =  trace (prestim + 1:prestim + 20);

asymptote = mean(y(16:20));
noise = 2*std(y(16:20));
maxAsym = asymptote + noise;
minAsym = asymptote - noise;

Ampli = max(y(1:4)) - asymptote;
maxAmpli = Ampli + noise;
minAmpli = Ampli - noise;




if numExp == 2
    F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata) + x(5);
    %% VARIABLES:
    % x(1) = amplitude of tau1
    % x(2) = 1 / tau1
    % x(3) = amplitude of tau2
    % x(4) = 1 / tau2
    % x(5) = asymptote
    if isempty(varargin)
        x0 = [Ampli 0.98 Ampli/2 5 asymptote] ;
        lb = [-inf -inf -inf -inf -inf];
        ub = [inf inf inf inf inf];
    else
        minTau = 1/varargin{1}; %% bigger
        maxTau = 1/varargin{2}; %% smaller

        %% minTau
        x0 = [max(y)*2 minTau max(y)*2 maxTau asymptote];
        
%         lb = [minAmpli maxTau/5 minAmpli minTau/5 minAsym];
%         ub = [maxAmpli minTau/5 maxAmpli minTau*5 maxAsym];
   %% median works well     
%         lb = [minAmpli 0 minAmpli 0 minAsym];
%         ub = [maxAmpli inf maxAmpli inf maxAsym];
%% decent mean, similar median than for the previous 
         lb = [minAmpli maxTau/5 minAmpli minTau/5 minAsym];
         ub = [maxAmpli minTau/5 maxAmpli inf maxAsym];
     
    
    end
    

    
    
else
    F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3);
    %% VARIABLES
    % x(1) = amplitude of tau1
    % x(2) = 1 / tau1
    % x(5) = asymptote
    if isempty(varargin)
        x0 = [Ampli 0.1 asymptote];
        lb = [-inf -inf -inf];
        ub = [inf inf inf];
        
%         x0 = [Ampli 0.1 asymptote];
%         lb = [minAmpli 0 minAsym];
%         ub = [maxAmpli inf maxAsym];
        
    else
        Tau = 1/varargin{1}; 
        x0 = [Ampli Tau asymptote];

        
        
        if Tau > 0
            %% 10 so far is  best
            lb = [minAmpli Tau/5 minAsym];
            ub = [maxAmpli Tau*5 maxAsym];
        else
            lb = [minAmpli 0 minAsym];
            ub = [maxAmpli inf maxAsym];
        
%             
%             lb = [minAmpli 0 minAsym];
%             ub = [maxAmpli inf maxAsym];
            
        end
        
%         


    end
    
end

%%

options = optimoptions(@lsqcurvefit,'Display','none','MaxIterations',500);%,'UseParallel', true);
[xunc, res]= lsqcurvefit(F, x0, t, y, lb , ub , options);
fitted =  F(xunc,t);

res = res/length(y);
if isempty(res)
    res = mean((y - fitted).^2);
end
%% normalize fitted and min sq

minFitted = min(fitted);
fittedN = fitted - minFitted;
maxFitted = max(fittedN);
fittedN = fittedN ./maxFitted;

yN = y - minFitted;
yN = yN ./maxFitted;
resN = mean((yN - fittedN).^2);




if fitted(1) > mean(fitted(2:3))
    ratio = (fitted(1) - mean(fitted(2:3)))/ fitted(1);
else
    ratio = 0;
end

% if ratio > 0.8
%     disp('bigOne')
% end


%% reorder and ratio: first the fast tau a its ratio, then the slow one

if numExp == 2
    [~,sl_loc] = max(xunc == min(xunc(2), xunc(4)));
    sl_tau = 1/xunc(sl_loc);
    sl_r =  xunc(sl_loc -1);
    x5 = xunc(5);

    if sl_loc > 3
        f_loc = 2;
    else
        f_loc = 4;
    end

    f_tau = 1/xunc(f_loc);
    f_r =  xunc(f_loc -1);

    y0 = fitted(1);
    f_r = f_r/ (y0 - x5);
    sl_r = sl_r/(y0 - x5);




    %% cases to refit with a single exp:
    %% 1) contribution of fast component is less than 10%
    %% 2) tauSl/tauFst < 2
    %% 3) any tau = 1000 seconds

    % if f_r >= 0.95 || f_r <= 0.05% || (sl_tau/f_tau) < 2 %|| sl_tau > 50 || f_tau > 5
    %     F = @(x,xdata)x(1)*exp(-x(2)*xdata)+ x(3);
    %     x0 = [1.4 0.5 -0.5] ;
    %     lb = [0 (1/20) min(y)/2];
    %     ub = [200 500 max(y)/2];
    %     [xunc, res] = lsqcurvefit(F, x0, t, y, lb , ub , options);
    %     
    %     
    %     fitted =  F(xunc,t);
    %     fm = mean(fitted(2:3));
    %     if fitted(1) > fm
    %         ratio = (fitted(1) - fm)/ fitted(1);
    %     else
    %         ratio = 0;
    %     end
    %     
    %     
    %     sl_tau = 1/xunc(2);
    %     sl_r = 1;
    %     f_tau = NaN;
    %     f_r =  NaN;
    %     x5 = xunc(3);
    %     
    %     if sl_tau > 50
    %        sl_tau = NaN;
    %        sl_r = NaN;
    %     end
    %     
    % end
     %% nov 2021 


    if f_r >= 0.95
        f_r = 1; 
        sl_r = 0;
        sl_tau = NaN;

    end
    if f_r <= 0.05
        sl_r = 1;
        f_r = 0;
        f_tau = NaN;

    end



    xunc(1,1) = f_r;
    xunc(1,2) = f_tau;
    xunc(1,3) = sl_r;
    xunc(1,4) = sl_tau;
    xunc(1,5) = x5;
    xunc(1,6) = res;
    xunc(1,7) = resN;
else
%     xunc(1,1) = f_r;
    xunc(1,2) = 1/xunc(1,2);
    xunc(1,4) = res;
    xunc(1,5) = resN;   
end


end

%a = 1, figure(1), plot(t, [averages.bin.Traces(:,a)]), hold on, plot(t, F([averages.bin.fit(a,:)],t);
