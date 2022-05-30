function rFatigue = rPfFagitue (Opts,pFMean,pFSen,b_v,nPar,yearV,SeaState)

    nYear   = numel(yearV);

    rFatigue = zeros(nPar*2, nYear);
    for ii = 1 : nYear
        yearSelect = yearV(ii); 

        [~,indexYear(ii)] = min(abs(Opts.yearsLifeExp - yearSelect )); 


        pFMeanIndex = pFMean (indexYear(ii), SeaState);

        % access data at the index position, each cell represents each
        % parameter
        pFSenIndex = cellfun(@(v) squeeze(v(:,indexYear(ii),SeaState)),pFSen);

        rFatigue(:,ii) = reshape(pFSenIndex.*b_v./pFMeanIndex, nPar*2, 1); % normalise , reshape to vector 

    end
end