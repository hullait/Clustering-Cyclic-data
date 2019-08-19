
%% colour sections of the test using labels

db = [161 4 211]./ 255;
brown = [220 119 26]./ 255;
Brown = [139 32 22]./ 255;
greeney = [40 112 8]./ 255;
purple= [140 62 214]./ 255;
citrine= [206 201 12]./ 255;

load('R93.mat')

load('L:\PassOff_Data\ClassificationAlgorithm\Lancaster\xwb_ETOPS_data.mat')
loc = isnan(data(:,4));
loc_p = find(loc==1);
y = data(:,2);

start = [1; loc_p];
endtime = [loc_p; length(y)];


        % sometimes the N1 speed contains negative values when engine is
        % not running, which we wil replace by zero
        for i=1:length(y)
            if y(i)<0
                y(i)=0;
            end
        end
        %label=class.labels;
        %start=class.start;
        %endtime=class.endtime;
        t=1:length(y);    
        h=figure();
        plot(t,y, 'k');
        hold on;
        for i=1:length(all_label)
            if all_label(i)==1
                index=start(i):endtime(i);
                %plot(t(index),y(index),'LineWidth', 1, 'Color', citrine);
                plot(t(index),y(index),'LineWidth', 1, 'Color', R(1+2,:));
                %txt1 = 'A';
                %text(index(1)+round(length(index)/2),max(y(index))+2,txt1)
            elseif all_label(i)==2
                index=start(i):endtime(i);
                %plot(t(index),y(index), 'b');
                plot(t(index),y(index),'LineWidth', 1, 'Color', R(2+2,:));
                %txt1 = 'B';
                %text(index(1)+round(length(index)/2),max(y(index))+2,txt1)
            elseif all_label(i)==3
                index=start(i):endtime(i);
                %plot(t(index),y(index), 'LineWidth', 1, 'Color', purple);
                plot(t(index),y(index),'LineWidth', 1, 'Color', R(3+2,:));
                %txt1 = 'C';
                %text(index(1)+round(length(index)/2),max(y(index))+2,txt1)
                %text(index(1)+round(length(index)/2),max(y(index)),txt1)
                
                %text(index(1)+round(length(index)/2),0,txt1)
            elseif all_label(i)==4
                index=start(i):endtime(i);
                %plot(t(index),y(index), 'c');
                plot(t(index),y(index),'LineWidth', 1, 'Color', R(4+2,:));
                %txt1 = 'R';
                %text(index(1)+round(length(index)/2),max(y(index))+2,txt1)
                
                %text(index(1)+round(length(index)/2),max(y(index)),txt1)
                %text(index(1)+round(length(index)/2),0,txt1)
            elseif all_label(i)==5
                index=start(i):endtime(i);
                %plot(t(index),y(index), 'LineWidth', 1, 'Color', greeney);
                plot(t(index),y(index),'LineWidth', 1, 'Color', R(5+2,:));
                %txt1 = 'P';
                %text(index(1)+round(length(index)/2),max(y(index)),txt1)
                %text(index(1)+round(length(index)/2),0,txt1)
              elseif all_label(i)==-1
                index=start(i):endtime(i);
                plot(t(index),y(index), 'r');
                %txt1 = 'U';
                %text(index(1)+round(length(index)/2),max(y(index)),txt1)
                %text(index(1)+round(length(index)/2),0,txt1)
            end
        end
        title('N1 speed plot of Cyclic Test')
        %title('N1 plot with epsilon=600')
        %title('N1 plot with labels using OPTICS')
        xlabel('time') % x-axis label
        ylabel('percentage speed') % y-axis label
        hold off;        
        %saveas(h, 'Cyclic_classification.png');
        %saveas(h, 'Cyclic_OPTICS_clusters_epsilon.png');
        saveas(h, 'Cyclic_clusters_epsilon1000.png');
