% calSen_KPIfree.m estimates the sensitivity using the eigenvectors of
% Fisher information 

% output 
%   - yjpdf: joint pdf of the random variables of y
%   - V_e  :  eigenvectors of fisher information matrix 
%   - D_e  :  eigenvalues


function [yjpdf,V_e,D_e] = calSen_KPIfree (y,xS,nPar,Ny,ListPar,parJ,isNorm,dummyVar)
     if nargin <= 7
        dummyVar = [];
     end

     yjpdf = cal_jpdf_hist (y,xS,Ny); % histogram approach

 
     if ~ isempty (dummyVar)
         % get rid of dummy variables 
         for ii = dummyVar(1) : dummyVar(end)
            for jj = 1 :4
                yjpdf.dp_y{ii,jj} = zeros(size(yjpdf.dp_y{ii,jj}));
            end
        
         end
     end

     % alternatively using CDF approach 
     % yjpdf = cal_pdf (y,xS,nPar);   % 1-D 
     % yjpdf= cal_jpdf_nD_2 (y,xS,nPar,Ny);  % n-D

     % compute fisher matrix (raw)
     Fraw = cal_jFisher (yjpdf,nPar);
     Fraw = Fraw(1:nPar*2,1:nPar*2);  % 3rd/4th parameters are not implemented  
     
     % re-parameterization and normalization  - transformation 
     Fn = parTran(Fraw, ListPar,parJ,isNorm) ;
     

     % eigen analysis
     [V_e,D_e] = eig(Fn); 
     
     lambda = diag(D_e);
     [EigSorter,EigIndex] = sort(lambda,'descend');
     V_e = V_e(:,EigIndex);
     lambda = lambda(EigIndex);
     D_e = diag(lambda);    