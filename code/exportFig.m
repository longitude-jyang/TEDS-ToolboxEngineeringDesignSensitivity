function exportFig(isExportFig,paperPath,figName)
    
    fileName = strcat(paperPath,'/',figName,'.png');
    if isExportFig == 1
        export_fig(fileName);
    end