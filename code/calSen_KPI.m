
function [pF,pFMean,pFSenC] = calSen_KPI (y,yExLevel,isPrctile,nQoI,xS)

   
    % size of random responses
    [Ns,N_QoI] = size(y); % Ns is number of samples
    
    pF = zeros(Ns,N_QoI);
    pFMean = zeros(1,N_QoI);
    pFSenC = cell(1,N_QoI);

    if nQoI <= N_QoI

        yExLevel = repmat (yExLevel,N_QoI);
 
    end

    % failure threshold 
    % could use either an absolute threshold, or a percentile (isPrctile == 1)
 
    for ii = 1 : N_QoI

        PeakThresholdFactor = yExLevel(ii) + (rand(1)-0.5); % add small noise to randomise the percentile every time the function is called

        if isPrctile == 1
              
            yPeakThreshold = prctile(y(:,ii) ,PeakThresholdFactor,'all'); % use percentile

        else
            yPeakThreshold = PeakThresholdFactor; % set absolute threshold for failure
        end

        yPeakThresholdGrid = repmat (yPeakThreshold,Ns,1);  % put the threshold into grid for comparison

        pF(:,ii) = y(:,ii) >= yPeakThresholdGrid;   % Indicator for failure        

        pFMean(ii) = mean(pF(:,ii), 1, 'omitnan'); % Unconditional probability of failure 
        
    
    % sensitivity
       
        [~,N_UPar]=size(xS.samp);
        pFSen = zeros(2,N_UPar);
        for kk = 1 : N_UPar
               pFSen (1,kk) = mean(pF(:,ii).*xS.senA(:,kk),1,'omitnan');  
               pFSen (2,kk) = mean(pF(:,ii).*xS.senB(:,kk),1,'omitnan'); 
        end 

        pFSenC{ii} = pFSen;  
    end
