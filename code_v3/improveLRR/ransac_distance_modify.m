function [o]=ransac_distance_modify(Xid,Params,TP)
%%
if TP.debug_class_study
    global r_positive r_negtive r_id;
end
nSample = size(Xid,2);

DEBUG = TP.debug_class;
numberall = TP.numberall;
kdtree = Params.kdtree;


for ri=1:nSample
    neighbours(:,ri) = kdtree_k_nearest_neighbors(kdtree,Params.Allpoints(Xid(ri),:),numberall)';
end
Params.Allpoints=Params.Allpoints';
mo = compute_similarity_matrix_a(Xid,neighbours,Params,TP);
o = modify_guiding_matrix(Xid, Params, TP, mo, r_positive, r_negtive, r_id);




