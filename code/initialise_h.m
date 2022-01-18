function [Opts, RandV] = initialise_h (CaseName,selpath_hfunction)
 

addpath (genpath(selpath_hfunction));

switch CaseName
    case { 'FWTtank'}
    %---------------------------
    % choose FWTtank 
    %---------------------------
    % (1) rho     = 1180
    % (2) rho_f   = 1025
    % (3) Ca      = 1
    % (4) Cd      = 1
    % (5) L_b     = 0.15
    % (6) r       = 0.045
    % (7) t       = 0.003
    % (8) mb      = 3
    % (9) E       = 3e9; 

    RandV.nVar = 9;
    RandV.varName=[{'\rho'},{'\rho_f'},{'C_a'},{'C_d'},{'L_b'},{'r'},{'t'},{'mb'},{'E'}]';
    
    RandV.vNominal = [1180 1025 1 1 0.15 4.5e-2 3e-3 3 3e9].';
    RandV.CoV = 1/1e2 * ones(RandV.nVar,1);

    Opts.funName  = 'design_FWTtank';
    Opts.distType = 'Normal';  
    
    
    case { 'Riser','riser'}
    %---------------------------
    % choose Riser
    %---------------------------   

    %     Opts.RanPar.Ca          = [1 1.5        1.5/5];
    %     Opts.RanPar.Cd          = [1 1.1        1.1/2];
    %     Opts.RanPar.rho         = [1 7.84e3     7.84e3/20];
    %     Opts.RanPar.E           = [1 2e11       2e11/20];
    %     Opts.RanPar.rho_f       = [0 1.025e3    0];
    %     Opts.RanPar.rho_oil     = [1 0.92e3     0.92e3/10];
    %     Opts.RanPar.T0          = [1 4905e3     4905e3/10];
    %     Opts.RanPar.Din         = [0 0.374      0];
    %     Opts.RanPar.Dout        = [0 0.406      0];
    %     Opts.RanPar.De          = [0 0.66       0];
    %     Opts.RanPar.snA         = [1 10^11.299  10^11.2990/2];
    %     Opts.RanPar.snB         = [1 3          3/2];

        % (1) Ca      = 1.5
        % (2) Cd      = 1.1
        % (3) rho     = 7.84e3
        % (4) E       = 2e11
        % (5) rho_oil = 0.92e3
        % (6) T0      = 4905e3


        RandV.nVar = 6;
        RandV.varName=[{'C_a'},{'C_d'},{'\rho'},{'E'},{'\rho_oil'},{'T_0'}]';

        RandV.vNominal = [1.5 1.1 7.84e3 2e11 0.92e3 4905e3].';
        RandV.CoV = [1/5 1/2 1/20 1/20 1/10 1/10].';

        Opts.funName  = 'design_Riser';
        Opts.distType = 'Normal';  
end