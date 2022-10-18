function [fitted, xunc]= gaussFit3CoCenter(trace, varargin)
%% varargin: fixed variables

t = (1:1:length(trace))';
[ln, colm] = size(trace);
if ln < colm
    y =  trace';
else
    y = trace;
end
%% a Gaussian plus a bckg = x(4)
F =  @(x,xdata)x(1)*exp(-((xdata-x(2))/x(3)).^2) + x(4);%*exp(-((xdata-x(2))/x(5)).^2);
%% free parameters
% x(1) = Amp 
% x(2) = loc peak
% x(3) = sigma
% x(4) = background

bck = mean(trace(1:2));

x0 = [max(trace) length(trace)/2 1.5 bck] ;
lb = [-inf 1 0.5 -inf];
ub = [inf length(trace) 10 inf];

if ~isempty(varargin)
    varargin = cell2mat(varargin);
    mskVariables = varargin ~= 0;
    varargin = varargin(mskVariables);
    x0(mskVariables) = varargin;
    lb(mskVariables) = varargin;
    ub(mskVariables) = varargin;
end






options = optimoptions(@lsqcurvefit,'Display','none','Diagnostics', 'off', 'MaxIterations',1000);%'UseParallel', true);
[xunc, ~]= lsqcurvefit(F, x0, t, y, lb , ub , options);
fitted =  F(xunc,t);


end