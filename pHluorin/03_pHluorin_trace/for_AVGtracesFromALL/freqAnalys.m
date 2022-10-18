function [f, P1] = freqAnalys(trace, prestim, lngth, varargin)
%X = (avgct19_a.Synch.bin(1).Traces(20:40,1))';
%trace = (avgct19_a.Synch.bin(4).meanTrace(20:40))';
%X = (avgct19_a.Synch.bin(1).meanTrace(40:59))';


trace = trace(prestim + 1:prestim + lngth);

Fs = 20;            % Sampling frequency Hz                  
%T = 1/Fs;             % Sampling period       
L = length(trace);             % Length of signal

%t = (0:L-1)*T;        % Time vector


Y = fft(trace);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (Fs*(0:(L/2))/L);

[a,b] = size(P1);
if a > b
    P1 = P1';
end
    
end


% figure(2), hold on, plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

