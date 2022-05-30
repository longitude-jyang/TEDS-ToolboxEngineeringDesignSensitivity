% callTEDS_riser estimates both KPI-free and KPI-based sensitivities for a
% design of marine riser. This is intended as a companion for the MSSP paper 
% Yang, J., Langley, R.S., Andrade, L., 2022. Digital twins for design in the presence of uncertainties. Mechanical Systems and Signal Processing 179, 109338. https://doi.org/10.1016/j.ymssp.2022.109338

% this is a seperate code from callTEDS, because instead of call the
% h_function in real time, the results from h_function are pre-saved to focus on sensitivity analysis
% the Monte Carlo data after h_function evaluation can be found in this
% dropbox folder: (https://www.dropbox.com/s/h5apdcgymz4yzfc/MR_RS2_FATIGUE_N5000_11-06-2020%2008-26.mat?dl=0)
% and h_function in this case uses the CHAOS code, which can be found from the repository 'https://github.com/longitude-jyang/hydro-suite'

% load data
%     load('MR_RS2_FATIGUE_N5000_11-06-2020 08-26.mat');
    PreCalc = R;
       
% ---------------------
% (0) options for Monte Carlo and Sensivity analysis 

    NsSelect = 5000;           % select number of samples 
    iUPar    = [1 2 3 4 6 7 11 12];  % index for uncertain parameters
    nPar     = numel(iUPar);
    
    
    % failure thresholds (KPI)
    fatigueLife = 80;   
    
    Ny = 40;       % number of bins for histogram of y pdf estimation 
    isNorm = 1;    % choose to normalise Fisher matrix (default == 1 for proportional normalization), 
                   %  2 for mean/std normalization,3 for std/std normalization 
   
    iState = 2; % seastate used 


% ---------------------
% (1) from data: distribution parameters
        [ParSen, b_v, ListPar,xS] = parDist(PreCalc,Opts,iUPar,NsSelect);   
        parJ = eye(nPar*2,nPar*2);  % identity if porportional normalization 

% (2) from data: quantities of interests (x3)    
        sigmaB    = squeeze(PreCalc.sigmaB(:,iState,:));  % stress
        Syy       = squeeze(PreCalc.Syy(:,iState,:));     % displacement
        Stheta    = squeeze(PreCalc.Stheta(:,iState,:));  % rotation 

% -------------------------------------------------------------------------
% (3)  sensitivity analysis using Fisher Information Matrix (KPI-free)  
% -------------------------------------------------------------------------
% with y, post process for FIM 
% first estimate the pdf of y and then form F matrix             
        
        % as the excitation is random, use peak r.m.s responses as QoIs        
        SyyPeak = max(Syy, [], 2);           
        SthetaPeak = max(Stheta, [], 2);
        sigmaBpeak = max(sigmaB, [], 2);
            
        y = [sigmaBpeak SyyPeak SthetaPeak ];

        % force the sensitivites to S-N parameters to zero(they are not involved in the calculation)
        % the effect on fisher is spurious 
        dummyVar = [7 8]; 

        disp(' KPI-free Sensitivity Analysis Starts: ...')
        tic;  
        
            [yjpdf,V_e,D_e] = calSen_KPIfree (y,xS,nPar,Ny,ListPar,parJ,isNorm,dummyVar);
            lambda = diag(D_e);
            V_e(:,lambda == 0) = 0; % if eigenvalue is zero (because of dummyVar), force eigenvector to zero

        
        elapseTime = floor(toc*100)/100; 
        disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))


% -------------------------------------------------------------------------
% (4)  sensitivity analysis for failure probability (KPI-based)  
% -------------------------------------------------------------------------   

     disp(' KPI-based Sensitivity Analysis Starts: ...')
     tic; 

        %  fatigue failure calls seperate code for KPI-based sensitivity,
        %  as there is random process involved 
        [Pf_fatigue,PfMean_fatigue,PfSen_fatigue] = calPfFatigue(PreCalc,Opts,iUPar,NsSelect,ParSen);

            pF_Fatigue.pF = Pf_fatigue;
            pF_Fatigue.pFm = PfMean_fatigue;
            pF_Fatigue.pFs = PfSen_fatigue;

         rFatigue = rPfFagitue (Opts,squeeze(pF_Fatigue.pFm),pF_Fatigue.pFs,b_v,nPar,fatigueLife,iState);    
         
     elapseTime = floor(toc*100)/100; 
     disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))        

