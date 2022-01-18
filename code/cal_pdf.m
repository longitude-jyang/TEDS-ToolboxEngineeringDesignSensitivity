function ypdf= cal_pdf (y,xS,nPar)

% (1) start from bending stress data for N random samples 

    [Ns,Ne]=size(y); % Ns is number of samples, Ne is the other dimension of y
    
    % obtain the bounds of y, note these could be vectors 
    yL=min(y); 
    yH=max(y);
    
% (2) compute cdf      
    Ny = 30;                           % define vector length of the cdf      <---- This is a variable
   
    dy=(yH-yL)/(Ny-1);                % step size vector 
    y_v=yL+[0:Ny-1].'*dy;             % y vector, at each position (Ny x Ne)
   
    
    % loop over the vector 
    P_y=zeros(Ny,Ne); % cdf vector of probability 
    dP_y=cell(nPar,4); % cdf vector of 1st derivative of probability
    
    
    for ii=1:Ny
        
        % loop over the samples, count how many of them are below threshold
        I_y=zeros(Ns,Ne); % indicator function, 2-D

        for jj=1:Ns        
                I_y(jj,y(jj,:)<=y_v(ii,:))=1;                
        end
        
        % average for expected probability 
        P_y(ii,:)=sum(I_y)/Ns;                       % for Py     
        
    
        for kk=1:nPar                                % for dPy
           
            dP_y{kk,1}(ii,:)=sum(I_y.*xS.senA(:,kk))/Ns;  
            dP_y{kk,2}(ii,:)=sum(I_y.*xS.senB(:,kk))/Ns;  
            dP_y{kk,3}(ii,:)=sum(I_y.*xS.senC(:,kk))/Ns;  
            dP_y{kk,4}(ii,:)=sum(I_y.*xS.senD(:,kk))/Ns;  
            
        end
        
    end

 % (3) compute pdf 
    p_y=diff(P_y)./diff(y_v);                       % length of Ny-1 x Ne
    p_y=[zeros(1,Ne);p_y];                          % length Ny x Ne
    dp_y=cell(nPar,4);
    
    for kk=1:nPar
           dp_y{kk,1}=[zeros(1,Ne); diff(dP_y{kk,1}) ./diff(y_v)];  
           dp_y{kk,2}=[zeros(1,Ne); diff(dP_y{kk,2}) ./diff(y_v)]; 
           dp_y{kk,3}=[zeros(1,Ne); diff(dP_y{kk,3}) ./diff(y_v)];  
           dp_y{kk,4}=[zeros(1,Ne); diff(dP_y{kk,4}) ./diff(y_v)];           
    end
    
% (4) assemble for output
    ypdf.p_y  = p_y;
    ypdf.dp_y = dp_y;
    ypdf.y_v  = y_v;
    ypdf.P_y  = P_y;