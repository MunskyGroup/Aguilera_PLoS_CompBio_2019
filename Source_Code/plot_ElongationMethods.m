function plot_ElongationMethods (folderName, namePlot, rsem_long, rsem_medium, rsem_short)

rsem_long (rsem_long>6)= 6;
rsem_medium(rsem_medium>6)= 6;
rsem_short(rsem_short>6) = 6;

maxCorrValue = max([max(max(rsem_long)),max(max(rsem_medium)),max(max(rsem_short))]);
maxCorrValue=6;

minCorrValue = min([min(min(rsem_long)),min(min(rsem_medium)),min(min(rsem_short))]);
minCorrValue =0;

if minCorrValue>0;minCorrValue=0;end
sizeWide=2.5;
sizeHeight=1;
labels_rep={'10','30','60','100'};

%% Plotting
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, sizeWide, sizeHeight];
hcb = imagesc(rsem_short);
colorbar
colormap(gray)
%colormap(flipud(gray));
set(gca,'YDir','normal') ;
caxis([minCorrValue maxCorrValue]);
xticks([1,2,3,4])
yticks([1,2,3,4])

matrix=rsem_short;
textStrings = num2str(matrix(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
for i =1:length(textStrings)
    if strcmp(textStrings{i},'6.00')==1
        textStrings{i} ='> 6.00';
    end
end
[x, y] = meshgrid(1:4);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
    'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(matrix(:) < midValue, 1, 3);  % Choose white or black for the
%   text color of the strings so
%   they can be easily seen over
%   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2),'FontSize',5);  % Change the text colors
xticklabels(labels_rep)
yticklabels({'1s','3s','10s','20s'})
xlabel('No. Spots','FontSize',12);
ylabel('Sampling Rate','FontSize',12);
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r600',[namePlot,'_short'])
movefile(horzcat([namePlot,'_short'], '.png'),horzcat(folderName),'f');


figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, sizeWide, sizeHeight];
hcb = imagesc(rsem_medium);
colorbar
colormap(gray)
%colormap(flipud(gray))
caxis([minCorrValue maxCorrValue]);
set(gca,'YDir','normal') ;
xticks([1,2,3,4]);
yticks([1,2,3,4]);
matrix=rsem_medium;
textStrings = num2str(matrix(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
for i =1:length(textStrings)
    if strcmp(textStrings{i},'6.00')==1
        textStrings{i} ='> 6.00';
    end
end
[x, y] = meshgrid(1:4);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
    'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(matrix(:) < midValue, 1, 3);  % Choose white or black for the
%   text color of the strings so
%   they can be easily seen over
%   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2),'FontSize',5);  % Change the text colors
xticklabels(labels_rep)
yticklabels({'1s','3s','10s','20s'})
xlabel('No. Spots','FontSize',12);
ylabel('Sampling Rate','FontSize',12);
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r600',[namePlot,'_medium'])
movefile(horzcat([namePlot,'_medium'], '.png'),horzcat(folderName),'f');



figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, sizeWide, sizeHeight];
hcb = imagesc(rsem_long);
colorbar
colormap(gray)
%colormap(flipud(gray));
caxis([minCorrValue maxCorrValue]);
set(gca,'YDir','normal')
xticks([1,2,3,4])
yticks([1,2,3,4])

matrix=rsem_long;
textStrings = num2str(matrix(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
for i =1:length(textStrings)
    if strcmp(textStrings{i},'6.00')==1
        textStrings{i} ='> 6.00';
    end
end
[x, y] = meshgrid(1:4);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
    'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(matrix(:) < midValue, 1, 3);  % Choose white or black for the
%   text color of the strings so
%   they can be easily seen over
%   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2),'FontSize',5);  % Change the text colors
xticklabels(labels_rep)
yticklabels({'1s','3s','10s','20s'})
xlabel('No. Spots','FontSize',12);
ylabel('Sampling Rate','FontSize',12);
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r600',[namePlot,'_long'])
movefile(horzcat([namePlot,'_long'], '.png'),horzcat(folderName),'f');


end