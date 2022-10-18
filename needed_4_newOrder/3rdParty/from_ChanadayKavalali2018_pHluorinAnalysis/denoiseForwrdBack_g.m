%% --- DENOISING OF THE FLUORESCENCE SIGNAL --- %
% Based on: Chung and Kennedy (1991); Reuel et al. (2001).
% It's a filtering on the time domain based on signal amplitudes.

% ^^^ NOTE: User can change following parameter to improve fit ^^^

function ROIclean = denoiseForwrdBack_g(ROI,M,P)
dnRounds=3;  % use 2 or 3;

%************************ LOOP in original
%nfile=1;
%in_filename='Data1.xlsx';
%out='Data1_';
%[status,sheets]=xlsfinfo(in_filename);
%NSheets=numel(sheets);

%************************ LOOP in original
%sheet=1;
%out_filename=[out,num2str(sheet)];
%Data=xlsread(in_filename,sheet);

%************************ PARAMETERS in original
%t = Data(:,3);
%F = Data(:,6:end);

%************************ BLEACHING and BACKGROUND correction in original
%ROI = F;


for round=1:dnRounds
    
N = [1 8];%[1 8];%1]; %8];% 16];% 32];% 64];  %Sampling window length (size of the forward and backward mean for the predictors)
K = length(N);      %  Number of forward and backward predictors
%M = 20; %M = 20;               Analysis window to compare the predictors
%P = 10; %P = 10;               Weighting factor

% --- Forward-Backward non-linear Algorithm

    for y = 1:size(ROI,2)
        testROI(:,1) = ROI(:,y);
        ltime=length(ROI);%ltime=length(t);
        I_avg_f = zeros(ltime,K);
        I_avg_b = zeros(ltime,K);
        
        for g = 1:ltime
            for k = 1:K
                % Average forward predictor
                window = N(k);
                if g == 1
                    I_avg_f(g,k) = testROI(1,1);
                elseif g - window - 1 < 0
                    I_avg_f(g,k) = sum(testROI(1:g-1,1))/g;
                else
                    epoint = g - window;
                    spoint = g - 1;
                    I_avg_f(g,k) = sum(testROI(epoint:spoint,1))/window;
                end
                % Average backward predictor
                if g == ltime
                    I_avg_b(g,k) = testROI(g,1);
                elseif g + window > ltime
                    sw = ltime - g;
                    I_avg_b(g,k) = sum(testROI(g+1:ltime,1))/sw;
                else
                    epoint = g + window;
                    spoint = g + 1;
                    I_avg_b(g,k) = sum(testROI(spoint:epoint,1))/window;
                end
            end
        end
        
        
        % Non-normalized forward and backward weights:
        f = zeros(ltime,K);
        b = zeros(ltime,K);
        for i = 1:ltime
            for k = 1:K
                Mstore_f = zeros(M,1);
                Mstore_b = zeros(M,1);
                Pi=1/(2*K); % Natali added this (see paper by Chung and Kennedy)
                for j = 0:M-1
                    t_f = i - j;
                    t_b = i + j;
                    if t_f < 1
                        Mstore_f(j+1,1) = (testROI(i,1) - I_avg_f(i,k))^2;
                    else
                        Mstore_f(j+1,1) = (testROI(t_f,1) - I_avg_f(t_f,k))^2;
                    end
                    % eqn. (4) and (5) in paper by Chung ang Kennedy:
                    if t_b > ltime
                        Mstore_b(j+1,1) = (testROI(i,1) - I_avg_b(i,k))^2;
                    else
                        Mstore_b(j+1,1) = (testROI(t_b,1) - I_avg_b(t_b,k))^2;
                    end
                end
                f(i,k) = Pi*(sum(Mstore_f)^(-P));
                b(i,k) = Pi*(sum(Mstore_b)^(-P));
            end
        end
        
        % Vector of normalization factors for the weights:
        C = zeros(ltime,1);
        for i = 1:ltime
            Kstore = zeros(K,1);
            for k = 1:K
                Kstore(k,1) = f(i,k) + b(i,k);
            end
            C(i,1) = 1/sum(Kstore);
        end
        
        % Putting parameters together and solving for intensities:
        ROIclean = zeros(ltime,1);
        for i = 1:ltime
            TempSum = zeros(K,1);
            for k = 1:K
                TempSum(k,1) = f(i,k)*C(i,1)*I_avg_f(i,k) + b(i,k)*C(i,1)*I_avg_b(i,k);
            % summatory over K of eqn. (2) in paper by Chung and Kennedy
            end
            ROIclean(i,1) = sum(TempSum);
        end
        
        ROI(2:ltime-1,y) = ROIclean(2:ltime-1,1);
        
    end

ROI(ROI==Inf)=500;
ROI(isnan(ROI))=0;

end
ROIclean = ROI;
%save(out_filename);
end
