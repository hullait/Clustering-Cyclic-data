
%% download drtoolbox

addpath('.\drtoolbox')

%% The Dynamic Time warping distances for the 281 manoeuvres in the Cyclic tests

dtw_distance = csvread('Cyclic_dtw_distances.csv');

%% K-nearest neighbour to choose epsilon for DBSCAN

[mIdx,mD] = knnsearch(dtw_distance,dtw_distance,'K',10,'Distance','euclidean');

    [outDTW, idDTW] = sort(mean(mD,2), 'descend');
    h = figure();
    plot(outDTW)
    %hline(1000, 'r')
    title('Ordered Minimum DTW')
    xlabel('Manoeuvres') % x-axis label
    ylabel('Min DTW')
    
%% DBSCAN    

epsilon = 1000;
min_clust = 10;
[C, ptsC, centres] = dbscan(dtw_distance, epsilon, min_clust);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Projections using drtoolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Projections - Laplacian 
%{
[mappedX, mapping] = compute_mapping(dtw_distance, 'Laplacian', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of Laplacian Eigenmaps'); drawnow

%% Projection - Isomap

[mappedX, mapping] = compute_mapping(dtw_distance, 'Isomap', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of Isomap'); drawnow

%% Projection - KernelPCA

[mappedX, mapping] = compute_mapping(dtw_distance, 'KernelPCA', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of KernelPCA'); drawnow

%% Projection - MDS

[mappedX, mapping] = compute_mapping(dtw_distance, 'MDS', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of MDS'); drawnow
%}

%% Projection - tSNE

ptsC(find(ptsC==0))=-1;
[mappedX, mapping] = compute_mapping(dtw_distance, 'tSNE', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of tSNE for Cyclic Data'); drawnow

%{
%% Projection - Sammon

[mappedX, mapping] = compute_mapping(dtw_distance, 'Sammon', 2);

%[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, gscatter(mappedX(:,1), mappedX(:,2), ptsC); title('Result of Sammon'); drawnow

%}

%csvwrite('tSNE_mapping_cyclic.csv',mappedX)
%csvwrite('cyclic_clusters', ptsC);