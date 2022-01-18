% cal_h is the call function to call the blackbox h-function and evaluate
% the quantity of interest (QoI), for all random samples 
% the size of y: [number of samples, number of QoIs]
% example of QoI: acceleration, stress, sound pressure etc 


% 22/04/2021 @ Franklin Court, Cambridge  [J Yang] --> initial trial 

function h_Results = cal_h (xS, Opts)

    % assign options
    nS = Opts.nSampMC;
    hfunction = str2func(Opts.funName);
    
    % do a run to check size of output y
    y0 = hfunction(xS.samp(1,:));
    Ne = length(y0);
    
    % loop over all samples 
    y = zeros(nS,Ne);
    
    waith = waitbar(0,'Please wait...');
    
    for ii=1:nS
        
        if strcmp(hfunction,'trail')                    % dummy function 
            y (ii,:)= sum(xS.samp(ii,:));
        end
        
        y(ii,:) = hfunction(xS.samp(ii,:));
        
        waitbar(ii / nS); 
    
    end
    
    close(waith)
    
    h_Results.y = y;