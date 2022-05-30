% parTran conducts re-parmaterization and normalizaiton, using matrix
% transformations 
% 23/09/2020 @ Franklin Court, Cambridge  [J Yang] 


function  Fn = parTran(Fraw, ListPar,parJ,isNorm) 
   if nargin <=3
       isNorm = 1;
   end


   % reparameterization to get back to mean and std, using Jacobian matrix 
   % in case of Gaussian, parJ is identity matrix 
    Fn = parJ.'*Fraw*parJ;  
    
   % normalization 
   if isNorm == 1   % proportional normalization 
       
        parA = ListPar (:,5); % para A 
        parB  = ListPar (:,6); % para B

        b_v = diag([parA;parB]); % b vector --> diagonal matrix 

        Fn = b_v*Fn*b_v; % normalization 
   
   elseif isNorm ==2   % mean/std normalization
       
        parMean = ListPar (:,2); % mean 
        parStd  = parMean.*ListPar (:,3); % std

        b_v = diag([parMean;parStd]); % b vector --> diagonal matrix 

        Fn = b_v*Fn*b_v; % normalization 
       
   elseif isNorm ==3   % std normalization 
       
       parMean = ListPar (:,2); % mean 
       parStd  = parMean.*ListPar (:,3); % std
        
       b_v = diag([parStd;parStd]); % b vector --> diagonal matrix 
    
       Fn = b_v*Fn*b_v; % normalization        
       
   end

