
%% download drtoolbox

addpath('.\drtoolbox')
addpath('.\drtoolbox\techniques')

%% The Dynamic Time warping distances for the 981 manoeuvres in the Pass-Off tests

% DTW distance matrix
dtw_distance = csvread('PassOff_dtw_distances.csv');

%% True labels of manoeuvres in Pass-Off tests

label  = csvread('PassOff_true_labels.csv');
number_man = length(label); % total number of maneouvress


%% Choose epsilon using k-nn search

[mIdx,mD] = knnsearch(dtw_distance,dtw_distance,'K',10,'Distance','euclidean');

[outDTW, idDTW] = sort(mean(mD,2), 'descend');
h = figure();
plot(outDTW)
hline(4000, 'r')
title('Ordered Minimum DTW')
xlabel('Manoeuvres') % x-axis label
ylabel('Min DTW')
%saveas(h, fullfile('L:\PassOff_Data\ClassificationAlgorithm\ETOPS\PassOff_knn_elbow.png'));

%% DBSCAN

epsilon = 4000;
min_clust = 10;
[C, ptsC, centres] = dbscan(dtw_distance, epsilon, min_clust);

num_groups = max(ptsC);
%% classification accuracy

% calculate precision and recall
true_label = label;
true_label(find(true_label==-1))=0;

number_clust = num_groups+1;
number_classes = 9;
num_in = ones(1,number_clust); % number of terms in each cluster
number_each_class = ones(number_classes, number_clust); % matrix of number of each class in each cluster
for jj= 1:number_clust
    clust1 = find(ptsC==(jj-1));
    num_in(jj) = length(clust1);
    cluster1 = ptsC(clust1);
    for ww=1:number_classes
        number_each_class(ww,jj) = length( intersect(clust1,find(true_label==(ww-1))) );
    end
end

number_true_classes_vec = sum(number_each_class, 2); % the true number of each class
each_tp = ones(1,number_clust);
each_fp = ones(1,number_clust);
for jj=1:number_classes
    if num_in(jj)==0
        return
    else
        each_tp(jj) = nchoosek(num_in(jj),2);
        each_fp(jj) = nchoosek(num_in(jj),2) - each_tp(jj);
    end
end

total = number_man*(number_man-1)/2;
true_pos = sum(each_tp);
false_pos = sum(each_fp);


% samples mislabelled as Unknown and one case of mislabelled RP manoeuvre,
%there are 17 manoeuvres we has labelled as Unknown that were put into their own cluster
false_neg = sum(number_each_class(2:end,1))+1+17;
true_neg = total -true_pos-false_pos-false_neg;

precision = true_pos/(true_pos+false_pos);
recall = true_pos/(true_pos+false_neg);

rand_term = (true_pos+true_neg)/total;

beta = 1; % weighting to penalise false negatives more than false positives
F_meas = (beta+1)*(precision*recall)/(beta^2*precision+recall);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Projections using drtoolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Projections - Laplacian 

% In Laplacian and Isomap the mapped data is smaller tha the number of samples

%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'Laplacian', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of Laplacian Eigenmaps'); drawnow
%}
%% Projection - Isomap
%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'Isomap',2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), label); title('Result of Isomap'); drawnow
%}
%% Projection - KernelPCA
%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'KernelPCA', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), label); title('Result of KernelPCA'); drawnow
%}
%% Projection - MDS
%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'MDS', 2);

figure, gscatter(mappedX(:,1), mappedX(:,2), label); title('Result of MDS'); drawnow % using true labels
%figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of MDS'); drawnow % using cluster labels
%}
%% Projection - tSNE - cluster labels
red = [255,0,0]./ 255;
load('R93.mat')

[mappedX, mapping] = compute_mapping(dtw_distance, 'tSNE', 2);

ptsC(find(ptsC==0))=-1;
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC, [red; R(1:5:90,:)], '..........', 10*ones(1,10)); title('Result of tSNE on Trent 1000 engines using cluster labels'); drawnow % using cluster labels

%tSNE_Trent1000_tests_cluster_labels.png
%% Projection - tSNE - true labels

word_label = string(true_label);

word_label(find(true_label==0))='U ';
word_label(find(true_label==1))='A ';
word_label(find(true_label==2))='B ';
word_label(find(true_label==3))='C ';
word_label(find(true_label==4))='R ';
word_label(find(true_label==5))='P ';
word_label(find(true_label==6))='RP';
word_label(find(true_label==7))='F ';
word_label(find(true_label==8))='V ';

red = [255,0,0]./ 255;
db = [161 4 211]./ 255;
blue = [0,0,255]./ 255;
purple = [128,0,128]./ 255;
cyan = [0,255,255] ./ 255;
brown = [220 119 26]./ 255;
citrine= [206 201 12]./ 255;
green = [0 255 0]./ 255;
magneta = [255 0 255]./ 255;

R_col = [db; blue; green; cyan; brown; magneta; purple; red; citrine];

%[mappedX, mapping] = compute_mapping(dtw_distance, 'tSNE', 2);

figure, 
h = gscatter(mappedX(:,1), mappedX(:,2), word_label, R_col, '.........', 9*ones(1,9)); 
title('Result of tSNE on Trent 1000 engines with true labels'); 
legend(h([8,1,2,7,4,5,9,3,6]), {'U ', 'A ', 'B ', 'C ', 'R ', 'P ', 'RP ', 'F ', 'V '});

%tSNE_Trent1000_tests_true_labels.png

%% Projection - Sammon
%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'Sammon', 2);

figure, gscatter(mappedX(:,1), mappedX(:,2), label); title('Result of Sammon'); drawnow
%figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of Sammon'); drawnow % using cluster labels
%}

