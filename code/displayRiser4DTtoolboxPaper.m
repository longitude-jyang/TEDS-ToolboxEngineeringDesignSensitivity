% this is the code to produce plots for DT toolbox paper 
% 26/10/2021 @ Franklin Court, Cambridge  [J Yang]

isExportFig = 0;
caseFolder = 'B';
figPath = 'Papers\DT sensitivity toolbox\figs';

varName = {'C_a', 'C_d', '\rho', 'E', '\rho_o', 'T_0','\alpha','\delta'};

%  (1)  
%  plot Fisher eigenvalues 
     
     fig1= figure;
     bar([1:nPar*2].',lambda)
     xlabel('Index of FIM EigValue')
     set(gca,'FontSize',14)
             xlim ([0.5 16.5])
          
    figName ='FisherEigValues';
    figuresize(24, 12, 'centimeters');
    movegui(fig1, [50 40])
    set(gcf, 'Color', 'w');
    
    exportFig(isExportFig,strcat(figPath,'\',caseFolder),figName);
    
% (2) 
% plot eigenvectors
     fig1 = figure;
     for ii = 1:6
         subplot(2,3,ii)
         b=bar([1:nPar*2]',V_e(:,ii));      
         set(gca,'FontSize',14)
         set(gca,'xtick',[round(nPar/2) nPar+round(nPar/2)],'xticklabel',[{'Mean'},{'Std Dev'}]);

         xtips = b.XData;
         ytips = b.YData;
         ytips = ytips.*double(ytips>0);
         labels = [varName varName];
         text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')

         ylim ([-1 1])
         title(['No.', num2str(ii) ,' FIM EigVector ','[\lambda_',num2str(ii),'=',num2str(round(lambda(ii)),'% 1.1e'),']'])
     end
     
     
    figName ='FisherEigvectors';
    figuresize(36, 18, 'centimeters');
    movegui(fig1, [50 40])
    set(gcf, 'Color', 'w');
    
    exportFig(isExportFig,strcat(figPath,'\',caseFolder),figName);
    

% (3)
% plot Pf sensitivity
  
    [~,indexYear] = min(abs(Opts.yearsLifeExp - fatigueLife ));
    pF1 = round(pF_Fatigue.pFm(1,indexYear,iState)*1e2)/1e2;

    % disp failure sensitivities
    titleStr = [{['Sensitivity-FatigueFailure ','[Pf = ',num2str(pF1),']']}];
    
        fig1 = figure;

        subplot(211)
        b = bar([1:nPar*2],rFatigue);         
        
        ttl = title(titleStr);
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
        ttl.HorizontalAlignment = 'left';  

        set(gca,'xtick',[round(nPar/2) nPar+round(nPar/2)],'xticklabel',[{'Mean'},{'Std Dev'}]);
        xtips = b.XData;
        ytips = b.YData;
        ytips = ytips.*double(ytips>0);
        labels = [varName  varName];
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    
        ylim([min(b.YData)-2 max(b.YData)+5])
    
        ylabel('r [-]')

        set(gca,'FontSize',14)
    

% (4)
    subplot(212)
    % projection to fisher eigenvectors, resutls from fisher calculation
    
    titleStr = [{' Projection onto FIM EigVector '}];
    
     V_fisher = V_e;
   
         
         rc = rFatigue;        
         s = rc.' * V_fisher / norm(rc);% 

         bar ([1:nPar*2] , abs(s));
 
        ttl = title(titleStr);
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
        ttl.HorizontalAlignment = 'left';  
         
  
        ylim ([0 1])
        xlim ([0.5 16.5])
        set(gca,'FontSize',14)
   
     
    figName ='PfFatigue';
    figuresize(24, 18, 'centimeters');
    movegui(fig1, [50 40])
    set(gcf, 'Color', 'w');    
    
    exportFig(isExportFig,strcat(figPath,'\',caseFolder),figName);