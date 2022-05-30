function [ListPar,J] = parList(Opts,RandV,isNorm)

% 05/05/2022 @ Franklin Court, Cambridge  [J Yang] --> add Gamma distribution 

% % assign options
    dist=Opts.distType;

    nVar = RandV.nVar;
    index = [1:nVar]';

    vNominal = RandV.vNominal;
    CoV = RandV.CoV;

    J = eye(nVar*2,nVar*2);  % identity if porportional normalization 

    switch dist
        case {'Normal'}
            distType = ones(nVar,1);  % No.1
            parA = vNominal;
            parB = vNominal.*CoV; 
            parC = NaN(nVar,1);
            parD = NaN(nVar,1);

            

        case {'Lognormal'} 
            distType = ones(nVar,1)*2;  % No.2 

            ParMean = vNominal; 
            ParStd = vNominal.*CoV; 

            mu = log(ParMean.^2./(sqrt(ParMean.^2+ParStd.^2)));
            sigma = sqrt(log(1+ParStd.^2./ParMean.^2));

            parA = mu;
            parB = sigma; 
            parC = NaN(nVar,1);
            parD = NaN(nVar,1);

            if isNorm == 2 || isNorm ==3 % if mean/std normalization, need to do re-parametrization 
                J = getJ (ParMean,ParStd,nVar,dist); 
            end

        case {'Gamma'}
            distType = ones(nVar,1)*3;  % No.3

            ParMean = vNominal; 
            ParStd = vNominal.*CoV; 

            k      = (ParMean./ParStd).^2 ; % shape parameter
            theta  = ParStd.^2./ParMean;  % scale parameter

            parA = k;
            parB = theta; 
            parC = NaN(nVar,1);
            parD = NaN(nVar,1);    

    end

    
    ListPar = [ index vNominal CoV distType parA parB parC parD];
    
%     for ii = 1:nVar
% 
%         ListPar(ii,:) = [ii vNominal(ii) CoV(ii) distType(ii) parA(ii) parB(ii) parC(ii) parD(ii)];
%     end


end

function J = getJ (m,s,nPar,dist)


    switch dist
        
         case {'Lognormal'}
    
            mu_m = 2./m - m./(m.^2 + s.^2);

            mu_s = - s./(m.^2 + s.^2);

            sigma_m = (log(1 + s.^2./m.^2)).^(-1/2).*(-s.^2./m)./(m.^2 + s.^2) ;

            sigma_s = (log(1 + s.^2./m.^2)).^(-1/2).*s./(m.^2 + s.^2) ;
    end


            J = zeros(nPar*2,nPar*2);

            J(1:nPar,1:nPar) = diag(mu_m);
            J(nPar+1:nPar*2,nPar+1:nPar*2) = diag(sigma_s);

            J(1:nPar,nPar+1:nPar*2) = diag(mu_s);
            J(nPar+1:nPar*2,1:nPar) = diag(sigma_m);
end