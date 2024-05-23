function AMBCAHeatmap(C_,savenclose_,figdir_)

C_norm = log2(C_); C_norm(C_norm < 0) = 0;
C_norm(logical(eye(size(C_norm)))) = 0;
C_norm = -(min(C_norm(:)) - C_norm) / (max(C_norm(:)) - min(C_norm(:))); 
figure('Units','inches','Position',[0 0 10 10]); 
imagesc(C_norm); colormap('hot'); clrbr = colorbar; axis square;
set(clrbr,'YTick',0:0.25:1);
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[],'FontName','Times','FontSize',20);

if savenclose_
    print([figdir_ filesep 'AMBCA_Heatmap'],'-dtiffn','-r300'); close;
end
end