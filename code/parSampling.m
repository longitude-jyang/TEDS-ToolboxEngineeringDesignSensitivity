function [xS,ListPar,ParSen] = parSampling(ListPar, nPar,nS)

% Sample of parameter list: (in columns)
% Name, NominalValue, Random, Distribution, A para, B para, C para, D para 
% 'E', 1e9, 1,'Normal',  1e9, 1e9/20                                        % <---- random example
% 'rho', 1, 0, [],                                                          % <---- deterministic example

% In this code, all parameters are considered to be random, if they are
% specified as design variables (determinstic), they will be given a delta
% approximated variance from uniform distribution                           % <---- issue with uniform: not continuous

% Samples are computed using the random function, which requires the name
% of the distribution considered and its corresponding parameters

xS.samp=zeros(nS,nPar);
xS.senA=zeros(nS,nPar);
xS.senB=zeros(nS,nPar);
xS.senC=zeros(nS,nPar);
xS.senD=zeros(nS,nPar);

ParSen = cell(nPar,2);

    for ii=1:nPar
        Par=ListPar(ii,:);
        
        if  Par(4)==1     
            dist='Normal';
            samp= random(dist, Par(5),Par(6), [1 nS]);   
            
        elseif Par(4)==2 
            dist='Lognormal';
            
            samp= random(dist, Par(5),Par(6),  [1 nS]);   

        elseif Par(4)==3 
            dist='Gamma';
            
            samp= random(dist, Par(5),Par(6),  [1 nS]); % shape, scale   
            
        else
            
            dist='Normal';
            mu=Par(2);
            st=Par(2)/1e4; 
            samp= random(dist, mu,st, [1 nS]);    
            ListPar(ii,5)=mu;                       % update with delta bound 
            ListPar(ii,6)=st;
            
            Par(5) = mu;
            Par(6) = st; 

        end


        % Sensitivity parameters (i.e. (dp/db)/p ) depend on the distribution 
        switch dist
            case {'normal', 'Normal'}
                mu = Par(5);
                st = Par(6);

                % Analytical differentiation for the normal distribution
                senA = (samp - mu)/st^2;
                senB= -1/st + (samp - mu).^2/st^3;    
                senC = NaN(1,nS);
                senD = NaN(1,nS);
                
            case {'Lognormal', 'LogNormal'}

                mu = Par(5);
                st = Par(6);
                
                
                % Analytical differentiation for the normal distribution
                senA = (log(samp) - mu)/st^2;
                senB= -1/st + (log(samp) - mu).^2/st^3;    
                senC = NaN(1,nS);
                senD = NaN(1,nS);

           case {'gamma', 'Gamma'}
                k = Par(5);
                theta = Par(6);

                % Analytical differentiation for the normal distribution
                senA = log(samp/theta) - psi(k); % psi is digamma function 
                senB = -k/theta + samp/theta^2 ;    
                senC = NaN(1,nS);
                senD = NaN(1,nS);


        end

        xS.samp(:,ii) = samp.';
        xS.senA(:,ii) = senA.';
        xS.senB(:,ii) = senB.';
        xS.senC(:,ii) = senC.';
        xS.senD(:,ii) = senD.';

        ParSen{ii,1} = senA.';
        ParSen{ii,2} = senB.';

    end
