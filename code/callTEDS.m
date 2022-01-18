% callTEDS.m calls the TEDS, toolbox for engineering design, and wraps
% around any blackbox functions/models/digital twins, to estimate the
% sensitivities to uncertainties

% 03/07/2020 @ Franklin Court, Cambridge  [J Yang] --> start with scalar parameter case
% 10/08/2020 @ Franklin Court, Cambridge  [J Yang] --> add option whether
% to normalise Fisher matrix with isNorm variable
% 22/04/2021 @ Franklin Court, Cambridge  [J Yang] --> in TEDS form 
% 04/05/2021 @ Franklin Court, Cambridge  [J Yang] --> add case for Riser 


% -------------------------------------------------------------------------
% (0)    
% -------------------------------------------------------------------------
% option section. in cases the toolbox is used for pre-saved results, can
% jump to section 3&4 for sensitivity analysis, with the adaption of the
% options from this section 
   
   % --------------------------------------------
   %  initialise the black box function h 

   selpath_hfunction = uigetdir([],'Select folder for h-functions:'); % select the directory for the h functions 


   % example cases, user defined 
   % the h-function name in the examples is defined as 'design_CaseName'
   CaseName = 'FWTtank' ;
   [Opts, RandV] = initialise_h (CaseName,selpath_hfunction);
    
   
    % --------------------------------------------
    % options for Monte Carlo and Sensivity analysis 

    Opts.nSampMC = 2000;    % number of samples for MC 
    Ny = 30;              % length y vector for cdf and pdf estimation 
    isNorm = 1;           % choose to normalise Fisher matrix (default == 1 for proportional normalization), 
                          %  2 for mean/std normalization,3 for std/std normalization 
    
    % options for KPI-based metric 
    nQoI  = 1;            % number of quantity of interest 
    yExLevel = 80;        % one threshold for each QoI 
    isPrctile = 1;        % is the threshold absolute or percentile?

% -------------------------------------------------------------------------
% (1)    
% -------------------------------------------------------------------------
% prepare random samples and store the distribution parameters in a matrix 'ListPar'


    [ListPar,parJ] = parList(Opts,RandV,isNorm);
    [nPar, ~] = size(ListPar);                                  % get size
    nS = Opts.nSampMC;                                          % No. of samples

    [xS,ListPar] = parSampling (ListPar, nPar,nS);   
   
% -------------------------------------------------------------------------
% (2)    
% -------------------------------------------------------------------------
% with the generated random samples, evaluate the blackbox function h 

     disp(' Monte Carlo Analysis Starts: ...')
     tic;  
     
        h_Results = cal_h (xS, Opts);
        y = h_Results.y;    
        [xS, y] = parFilter (xS, y);
        
     elapseTime = floor(toc*100)/100; 
     disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))
     
% -------------------------------------------------------------------------
% (3)  sensitivity analysis using Fisher Information Matrix (KPI-free)  
% -------------------------------------------------------------------------
% with y, post process for FIM 
% first estimate the pdf of y and then form F matrix     
     

     disp(' KPI-free Sensitivity Analysis Starts: ...')
     tic;  

     [yjpdf,V_e,D_e] = calSen_KPIfree (y,xS,nPar,Ny,ListPar,parJ,isNorm);
     
     elapseTime = floor(toc*100)/100; 
     disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))
     
    % ---------------------------------------
    % show sensitivity results 

    displaySen_KPIfree (D_e,V_e,nPar,RandV)
         
         
% -------------------------------------------------------------------------
% (4)  sensitivity analysis for failure probability (KPI-based)  
% -------------------------------------------------------------------------   

     disp(' KPI-based Sensitivity Analysis Starts: ...')
     tic;  
         [pF,pFMean,pFSen] = calSen_KPI (y,yExLevel,isPrctile,nQoI,xS);
         
     elapseTime = floor(toc*100)/100; 
     disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))
