function y = design_FWTtank (xS)

% this is called by cal_h for each sample 
% adjust the QoI outputs here 



     % get parameter values for the input, including the random ones 
     ModPar = getModPar_FWTtank (xS) ;   
    

          % call maincode.m to evaluate the h-function 
          %---------------------------------------------------------------------    
            Out = maincode (ModPar);    
          %---------------------------------------------------------------------

 
     % pick up response 
     response = Out.response;
     BM       = Out.BM;
     wave     = Out.wave;
     mode     = Out.mode;
     
     om_range = wave.F_f;
     [~,iom]  = min(abs(om_range-sqrt(mode.om(2))));  
     
     
    
% %      % interested in displacement response 
% %       z  = abs(response.frf); 
%      
%        % interested in stress response (rigid body)
%        z = abs(BM.SigmaIR);
%        % peak QoIs 
%        y = max(z (:,iom));


       % interested in natural frequencies 
       z = sqrt(mode.om(2));  % om was in omega^2  [rad/s]

       y = real(z(:).');   % for floating stucture, natural frequency can be negative for unstable structures  
       


     
     
 
 
 
 
    
    
    