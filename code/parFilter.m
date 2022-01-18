function [xS, y] = parFilter (xS, y)
   

    % remove NaN or infinites 
    index = prod(isfinite(y),2);
    
    index = logical(index);
    
    y = y(index, :);
    


    samp = xS.samp;
    senA = xS.senA;
    senB = xS.senB;
    senC = xS.senC;
    senD = xS.senD;
    
    xS.samp = samp(index,:);
    xS.senA = senA(index,:);
    xS.senB = senB(index,:);
    xS.senC = senC(index,:);
    xS.senD = senD(index,:);
    
    % remove outliers
   [y ,TF] = rmoutliers(y,'percentiles',[1 99]); 
    
  
    xS.samp = xS.samp(~TF,:);
    xS.senA = xS.senA(~TF,:);
    xS.senB = xS.senB(~TF,:);
    xS.senC = xS.senC(~TF,:);
    xS.senD = xS.senD(~TF,:);