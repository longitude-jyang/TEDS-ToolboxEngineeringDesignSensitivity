function [pF,pFMean,pFSen] = calPfyEx(yV,Opts,dispExLevel,NsSelect,ParSen,mode)
   
    %  displacement 
    yPeak = max(yV,[],2); % get peak, as a vector for all samples 
        

    % failure threshold 
        nThresholdFactor = numel(dispExLevel);
        PeakThresholdFactor = dispExLevel ; %+ (rand(1)-0.5) add small noise to randomise the percentile every time the function is called

        pF  = zeros(NsSelect,nThresholdFactor);

        for ii = 1 : nThresholdFactor

            if mode == 1
                yPeakThreshold = prctile(yPeak ,PeakThresholdFactor(ii)); % use percentile
            elseif mode ==2 
                yPeakThreshold = PeakThresholdFactor(ii); % use percentile

            end

            yPeakThresholdGrid = repmat (yPeakThreshold,NsSelect,1);  % put the threshold into grid for comparison

            pF(:,ii) = yPeak >= yPeakThresholdGrid;   % Indicator for failure        

        end
        pFMean = mean(pF, 1, 'omitnan');% Unconditional probability of failure 
        
        if iscolumn(pFMean)
            
            pFMean = pFMean.';  % if only one threshold, the squeeze gives column vector; put it as row vector
        end
    
    % sensitivity
        pFSen = cell(size(ParSen));
        [N_UPar,~]=size(ParSen);
        for kk = 1 : N_UPar
               pFSen {kk,1} = mean(pF.*ParSen{kk,1},1,'omitnan');  
               pFSen {kk,2} = mean(pF.*ParSen{kk,2},1,'omitnan'); 
        end 