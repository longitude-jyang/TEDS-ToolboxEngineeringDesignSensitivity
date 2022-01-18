
function y = design_Riser (xS)

     ModPar = getModPar_Riser (xS) ;   
    
          %---------------------------------------------------------------------    
            Out = maincode (ModPar);    
          %---------------------------------------------------------------------

 
 
    % pick up response 
     response = Out.response;
     BM       = Out.BM;
     wave     = Out.wave;
     mode     = Out.mode;
     

    % for random wave, r.m.s response computed for the riser and the output
    % is spatial dependent 

    rms_SmmEI   =  BM.rms_SmmEI   ;          % r.m.s bending moment 
    rms_SigmaEI =  BM.rms_SigmaEI ;          % r.m.s stress 

    rms_Syy     =  response.rms_Syy ;        % r.m.s displacement

 
    SyyPeak     = max(rms_Syy);           
    SigmaPeak   = max(rms_SigmaEI);
    
    % output QoIs 
    y = [SyyPeak SigmaPeak];
    