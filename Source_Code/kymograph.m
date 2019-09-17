function kymograph (s_out,nameplot,folderName)
s_out =s_out(1:500,:);
for i =1:size(s_out,1)
    for j =2:size(s_out,2)-10
        if s_out(i,j-1) ==0 && s_out(i,j)==1
   s_out(i,j:j+9)=1;
        end
    end
end

%% Plotting
figure('visible', 'off');
fig1= gcf;
 fig1.PaperUnits = 'inches';
 fig1.PaperPosition = [0, 0, 5.5, 6.1];
hold on
imshow(s_out)
axis tight
 h = gca;
h.Visible = 'On';
box on
set(gca,'linewidth',2)
set (gca ,'FontSize',12); 
set(gca, 'FontName', 'Arial');
xlabel('Codon Position','FontSize',14,'FontName', 'Arial');
ylabel('Time (sec)','FontName', 'Arial','FontSize',14);
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

end
