    
%%
tSNE_map  = csvread('tSNE_mapping_cyclic.csv');
clusters = csvread('cyclic_clusters.csv');

X = tSNE_map(:,1);
Y = tSNE_map(:,2);

N = length(X);

%%
%{
minX = min(X);
maxX = max(X);

minY=min(Y);
maxY=max(Y);

figure;
%colormap(customMap);
axes('position',[0 0 1 1]);
plot1 = scatter(X(1),Y(1),30,clusters(1),'.');
xlim([lonMin lonMax]);
ylim([latMin latMax]);
set(gca,'Color','none');
set(gca,'CLim',[0, 1E-4]);

for k = 2:length(lat) 
     plot1.XData = X(1:k); 
     plot1.YData = Y(1:k); 
     plot1.CData = clusters(1:k); 
     % pause 2/10 second: 
     pause(0.2)
end
%}
%%

% colours
%load('R93.mat')

minX = min(X)-2;
maxX = max(X)+2;

minY=min(Y)-2;
maxY=max(Y)+2;

cmap = hsv(5); %// define your colormap here
for i = 1:N
    figure(1)    
    xlim([minX maxX]);
    ylim([minY maxY]);
    %set(gca, 'Position', [0.23 0.2 0.18 0.65])
    %set(gcf,'Position',[60 502 1000 220])
    gscatter(X(1:i),Y(1:i), clusters(1:i), cmap);
    hold on
    %hold off
    F(i) = getframe(gcf) ;
    drawnow
end

% create the video writer with 1 fps
writerObj = VideoWriter('myVideo.avi');
writerObj.FrameRate = 7;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);