% cal_jFisher calculates the fisher matrix (raw)
% 29/09/2020 @ Franklin Court, Cambridge  [J Yang] --> removed
% normalization steps. normalization is now done with transformation and
% reparameterization 


function F = cal_jFisher (yjpdf,nPar)


    % for joint pdf, these are all Ny x Ny matrix
    p_y  = yjpdf.p_y;
    dp_y = yjpdf.dp_y;
    y_v  = yjpdf.y_v;

    % initialise vectors 
    [Ny,Ne]=size(y_v);
    
    F=zeros(nPar*4,nPar*4);
    
    % compute the grid size for Fisher intergration 
    dy_v=[zeros(1,Ne);diff(y_v)]; 
    
    
    % to allow different grid size (variable for bin size of jpdf), and
    % high deminsions, use cells and take products of multiple cells for
    % correct grid size --> this may become a problem with large matrix
    
        %  creat Ne dimensional grid for dy_v  
           dy_v_cell = num2cell (dy_v, 1)  ;  % each column goes to a cell 
           [dy_v_ndgrid{1:numel(dy_v_cell)}] = ndgrid(dy_v_cell{:});

        dynD = 1;    
        for ii = 1 :Ne
            dynD = dynD .*cell2mat(dy_v_ndgrid(ii)); % grid size matix, multiply by each dimension at a time
        end
 

    % loop over the parameters 
    for ii=1:nPar*4

        if ii<=nPar        
            
            dp_i=dp_y{ii,1};

        elseif ii<=nPar*2

            dp_i=dp_y{ii-nPar,2};
 
        elseif ii<=nPar*3

            dp_i=dp_y{ii-nPar*2,3};
            
        elseif ii<=nPar*4

            dp_i=dp_y{ii-nPar*3,4};
        end


        for jj=1:nPar*4
            if jj<=nPar

                dp_j=dp_y{jj,1};

            elseif jj<=nPar*2

                dp_j=dp_y{jj-nPar,2};

            elseif jj<=nPar*3

                dp_j=dp_y{jj-nPar*2,3};

            elseif jj<=nPar*4

                dp_j=dp_y{jj-nPar*3,4};
            end

            
            Fintegrand = dp_i.*dp_j./p_y.*dynD; 

            F(ii,jj) = sum(Fintegrand(isfinite(Fintegrand)),'all');  % non-normalised


        end
    end