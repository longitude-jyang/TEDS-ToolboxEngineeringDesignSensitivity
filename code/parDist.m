% parDist calculates the derivative term analytically for each sample 

function [ParSen, b_v, ListPar, xS] = parDist(PreCalc,Opts,iUPar,NsSelect)

% read distribution parameters
    RanPar = Opts.RanPar;
    UPar = struct2cell(RanPar);

    % read sample data
    parSet = PreCalc.parSet;

    N_UPar  = length(iUPar);

    ParSen = cell(N_UPar,2);
    b_v    = zeros(N_UPar,2);
    
    ListPar = NaN (N_UPar , 8);
    for kk = 1 : N_UPar

            Par = UPar{iUPar(kk)};
            mu  = Par(2);
            st  = Par(3);
            
            ListPar (kk, [5 6]) = [mu st];

            samp = parSet(1 : NsSelect,iUPar(kk));

            % Analytical differentiation for the normal distribution
            senA = (samp - mu)/st^2;
            senB = -1/st + (samp - mu).^2/st^3; 
            senC = zeros(size(senA));
            senD = zeros(size(senA));

            ParSen{kk,1} = senA;
            ParSen{kk,2} = senB;
            b_v(kk,:) =[mu st];
            
            
            xS.samp(:,kk) = samp; 
            xS.senA(:,kk) = senA;
            xS.senB(:,kk) = senB;
            xS.senC(:,kk) = senC;
            xS.senD(:,kk) = senD;
    end    
