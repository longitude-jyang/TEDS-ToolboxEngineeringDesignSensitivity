% cal_jpdf_hist  estimates the joint pdf using the method of histogram 
% for the purpose of Fisher information matrix, it estimates not only the
% jpdf (p(y)) but also its sensitivity (dp/db) 
% found this way to be much more quicker than diff of the jcdf 

% 02/12/2021 @ Franklin Court, Cambridge  [J Yang] 

function yjpdf = cal_jpdf_hist (y,xS,Ny)

    % test data
    % y = randn(50,2);

    [Ns,Ne] = size(y); % Ns is number of samples, Ne jpdf dimension 
      nbins = Ny - 1;      % number of bins for all dimensions
     
    edge_v  = num2cell(ones(1,Ne)*nbins); 

    % uses the histcn.m from matlab fileExchange, with slight modifications
    % fit for purpose for now, but still needs to improve it --------
    [epdf,epdf_dp, edges, h] = histcn(y,xS, edge_v{:});

    % convert the edges to matrix 
        y_v= cell2mat(edges.');


    % assemble for output
    yjpdf.p_y  = epdf;
    yjpdf.dp_y = epdf_dp;
    yjpdf.y_v  = y_v.';




