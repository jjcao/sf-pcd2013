%   Distribution code by JJ Cao Copyright 2013.
%
%   The Code is created based on the method described in the following paper 
%   [1] "Point Cloud Normal Estimation via Low-Rank Subspace Clustering", Jie Zhang, Junjie Cao, Xiuping Liu, Jun Wang,  Jian Liu, Xiquan Shi
%   (Shape Modeling International 2013), 2013. 
%  
%   The code and the algorithm are for non-comercial use only.

%% setup toolbox
clear;clc;close all;
addpath('toolbox/cvx')
addpath ('toolbox/jjcao_io')
addpath ('toolbox/jjcao_point')
addpath ('toolbox/jjcao_common')
addpath ('toolbox/jjcao_interact')
addpath ('toolbox/jjcao_plot')
addpath ('toolbox/kdtree')

addpath('improveLRR')
addpath('toolbox/zj_fitting')
cvx_setup

%% debug options
ADDNOISE = 0;

TP.debug_data = 0;
TP.debug_class = 0;
TP.debug_class_study = 1;
TP.debug_inclass_number = 1;
TP.debug_feature_study = 1;

%% algorithm options
TP.k_knn_feature = 70; % tunable arguments
TP.k_knn_normals = 30;
TP.k_knn_origin = 15;
TP.k_knn  = 120;
TP.Nneighbors = 100;

TP.sigma_threshold = 0.05;

TP.numberall = 30; % TP.numberall > TP.k_knn_normals;
TP.numberi =10;
TP.guide_rate = 0.4;
TP.tolerate_angle = 0.7; % no greater than pi/2

TP.lambda = 1;     %%ÔëÉù
TP.gama = 1;

%% read input && add noise && plot it && build kdtree
[P.pts] = read_noff('test.off');

% add noise
if ADDNOISE
    type = 'gaussian';%type = 'random';% type = 'gaussian';% type = 'salt & pepper';
    base = 'average_edge';%base = 'average_edge'% base = 'diagonal_line';
    p3 = 0.5;
    P.kdtree = kdtree_build(P.pts);
    P.pts = pcd_noise_point(P.pts, type, base, p3,P.kdtree);
    kdtree_delete(P.kdtree);
end

% build kdtreeÓÃÀ´ÕÒk½üÁÚ
P.kdtree = kdtree_build(P.pts);
TP.ever_length = compute_average_radius(P.pts,2,P.kdtree);

if TP.debug_data;
    figure('Name','Input'); set(gcf,'color','white');set(gcf,'Renderer','OpenGL');
    movegui('northeast');
    scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),30,'.','MarkerEdgeColor',GS.CLASS_COLOR5);  hold on;
    axis off;axis equal;
    view3d rot; % vidw3d zoom; % r for rot; z for zoom;
end

%% compute initial features (a ribbon)
[sigms normVectors errs]= compute_points_sigms_normals(P.pts, TP.k_knn_feature, P.kdtree);

TP.feature_threshold = feature_threshold_selection(sigms,TP);
P.init_feat = find(sigms > TP.feature_threshold);
%  the choose of fitting threshold and inner
TP.fitting_threshold = fitting_threshold_selection(errs,P.init_feat);
TP.inner_threshold = 2 * TP.fitting_threshold;
disp('initial features: ');

%% order of neighbor cluster, we wish compute edge_feature first
feature_sigms = sigms(P.init_feat);
mean_feature_sigms = mean(feature_sigms);
di_sigms = abs(feature_sigms - mean_feature_sigms);
[value_sigms, id_sigms] = sort(di_sigms);
id_feature = P.init_feat(id_sigms);
disp('order of cluster: ');


nFeature = length(id_feature);

%% initialize variable
nSample = size(P.pts,1);
labels = cell(1,nFeature);
feature = [];
TP.feature_indicator = zeros(1,nSample);
if TP.debug_class_study
    global r_positive r_negtive r_id;
    r_positive = zeros(nFeature,2*TP.k_knn);
    r_negtive = zeros(nFeature,2*TP.k_knn);
    r_id = zeros(nFeature,2*TP.k_knn);
    for j = 1:nFeature
        r_id(j,:) = kdtree_k_nearest_neighbors(P.kdtree,P.pts(id_feature(j),:),2*TP.k_knn);
        TP.feature_indicator(id_feature(j)) = j;
    end
end
P.normals = compute_normals(P.pts, TP.k_knn_normals, P.kdtree); % for construct a 6-dimension vector
Params.sigms = sigms;
Params.kdtree = P.kdtree;
Params.Allpoints = P.pts;
Params.normvectors = P.normals;
disp('initialize variable: ')

%% compute nomals 
for i = 1:nFeature;
    ii = id_feature(i);
    Params.current = P.pts(ii,:)';
    
    knn = kdtree_k_nearest_neighbors(P.kdtree,P.pts(ii,:),TP.k_knn)';
    knn_origin = kdtree_k_nearest_neighbors(P.kdtree,P.pts(ii,:),TP.k_knn_origin)';
    

    %% compute original point
    if sum(sigms(knn_origin))==0
        weight = 1/TP.k_knn_origin;
    else
        weight = sigms(knn_origin)/sum(sigms(knn_origin));
    end
    Params.origin = P.pts(knn_origin,:)'*weight;
    
   %%  Calculation subspace data via original point
    tX = (P.pts(knn,:) - repmat(Params.origin',TP.k_knn,1))';
    normtX = sum(tX.*tX);
    [zero zero_id] = min(normtX);
    tX(:,zero_id) = [];  normtX(zero_id) = [];
    tX1 = tX./repmat(normtX,3,1);
    knn(zero_id) = [];
    
   %%   classification
    tN = P.normals(knn,:)';
    X = [tX;tN];
    [norm_vector] = improveLRR_feature(X,knn,Params,TP);
    normVectors(ii,:) = norm_vector;
end
%%
kdtree_delete(P.kdtree);
write_noff('test_result.off',P.pts,  normVectors)

