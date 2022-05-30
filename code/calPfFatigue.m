
% calPfFatigue calculates the fatigue failure 
function [pF,pFMean,pFSen] = calPfFatigue(PreCalc,Opts,iUPar,NsSelect,ParSen)

    % bending stress 
    sigmaB      = PreCalc.sigmaB(1 : NsSelect,:,:);   
    sigmaBDot   = PreCalc.sigmaBDot(1 : NsSelect,:,:);
    
    [Ns,NT,Ne]  = size(sigmaB); % Ns is number of samples, Ne is number of elements along the riser, NT is number of sea states
    
    % Probability of occurrence of each sea storm (precomputed information).
    % The repmat is used for a later elementwise multiplication
    
    pSeaState = repmat(Opts.DetPar.jonswapPSS(1:Opts.DetPar.jonswapNSS)', [1 Ne]);
    
    % read S-N parameters 
    parSet = PreCalc.parSet;
    
    snARan = parSet(1 : NsSelect,iUPar(end-1));
    snBRan = parSet(1 : NsSelect,iUPar(end));
    
    
    % fatigue failure 

        % Design service time considered (in seconds) 
        tLife = Opts.yearsLifeExp*(365*24*3600);
        
        % Number of life times considered
        NL = length(tLife);
        
        % frequency for bending stress cycles
        freqCycle = (1/(2*pi))*(sigmaBDot./sigmaB);
        
        % initialise pF matrix 
        pF = zeros(Ns, NL, NT + 1 );
        
        for iS = 1 : Ns 
            
            % Damage per second with random parameters
            dmgSecRan = squeeze(freqCycle(iS,:,:)./snARan(iS).*...
                (sqrt(2)*sigmaB(iS,:,:)).^snBRan(iS).*gamma((snBRan(iS)+2)/2));
            
            % Thresholds loop
            for iL = 1 : NL                
                
                % Total damage recieved by the marine riser 
                % put results from individual sea states, first NT
                % rows, and sum of sea states, last row, in the same
                % matrix
                totDmg =[dmgSecRan ; ...
                    sum(dmgSecRan.*pSeaState, 1)]*tLife(iL);
                
                % Maximum damage along the marine riser length
                maxDmg = max(totDmg, [] ,2); % take max along columns 
                 
                % Total damage multiplied by the design fatigue failure
                dmgWithDff = maxDmg*Opts.failureDFF;
        
                % Probability of failure cases:
                % Deterministic material fatigue capacity equal to 1
                for iT = 1 : NT + 1 
                    if dmgWithDff (iT) > 1
                        pF(iS,iL,iT) = 1;
                    else
                        pF(iS,iL,iT) = 0;
                    end
                end
                
            end
            
        end
        
        % Unconditional probability of failure
        pFMean = mean(pF, 1, 'omitnan');
        
        
        
        % sensitivity
        pFSen = cell(size(ParSen));
        N_UPar  = length(iUPar);
        for kk = 1 : N_UPar
               pFSen {kk,1} = mean(pF.*ParSen{kk,1},1,'omitnan');  
               pFSen {kk,2} = mean(pF.*ParSen{kk,2},1,'omitnan'); 
        end   