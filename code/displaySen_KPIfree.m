function displaySen_KPIfree (D_e,V_e,nPar,RandV)

 format short e

 varName = RandV.varName;
 
 
 % plot Fisher eigenvalues 
 lambda = diag(D_e);
 
 figure
 bar([1:nPar*2].',lambda)
 xlabel('Index of Fisher Eigvalue')
 set(gca,'FontSize',14)
 
 % plot eigenvectors
 figure;
 for ii=1:4
     subplot(2,2,ii)
     b=bar([1:nPar*2]',V_e(:,ii));      
     ylim([-1 1])
     
     set(gca,'FontSize',14)
     set(gca,'xtick',[round(nPar/2) nPar+round(nPar/2)],'xticklabel',[{'Mean'},{'Std Dev'}]);

     xtips = b.XData;
     ytips = b.YData;
     ytips = ytips.*double(ytips>0);
     labels = [varName;varName];
     text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
    
     title(['No.', num2str(ii) ,' Fisher Eigevector ','[\lambda_',num2str(ii),'=',num2str(round(lambda(ii)),'% 1.1e'),']'])
 end
 
 
 

 
 