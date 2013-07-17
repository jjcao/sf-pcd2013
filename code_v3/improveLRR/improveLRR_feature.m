function [norm_vector label] = improveLRR_feature(X,Xid,Params,TP)

if TP.debug_class_study
    global r_positive r_negtive r_id;
end
DEBUG = TP.debug_class;
nSample = size(X,2);

%%   compute identity matrix
[o]=ransac_distance_modify(Xid,Params,TP);

%%
nGroup = 2;
[L,E] = admreglrr(X,TP.lambda,o,TP.gama);

for iRecursion = 1:4
    label = clu_ncut(L,nGroup);
    if DEBUG
        plot_claasification_points(Params.Allpoints, Xid, label);
    end
    tempX = Params.Allpoints(Xid(1:TP.Nneighbors),:)'; tempLabel = label(1:TP.Nneighbors);
    tempIds = cell(1,nGroup); tempXs = cell(1,nGroup);
    disting = ones(1,nGroup);
    for iGroup = 1:nGroup
        tempIds{iGroup} = find(tempLabel == iGroup);
        tempXs{iGroup} = tempX(:,tempIds{iGroup});
        if length(tempIds{iGroup})<3
            disting(iGroup) = 1;
        else
            [disting(iGroup)] = pca_fitting(tempXs{iGroup},TP);
        end
    end
    n_nonflat = length(find(disting==0));
    
    if (n_nonflat==0)||(iRecursion == 4)  
        if (n_nonflat==0)
           rlabel = zeros(1,TP.Nneighbors);
            for iGroup = 1:nGroup
                if length(tempIds{iGroup})<3
                    rlabel(tempIds{iGroup}) = 0;
                else
                    rlabel(tempIds{iGroup}) = iGroup;
                end
            end
            for ri = 1:TP.Nneighbors
                if rlabel(ri) == 0 || Params.sigms(Xid(ri)) < TP.feature_threshold
                   continue
                end
                tri = TP.feature_indicator(Xid(ri));
                if tri == 0
                   continue
                end
                for ci = 1:TP.Nneighbors
                    if rlabel(ci) == 0
                        continue
                    end
                    tci = find(r_id(tri,:) == Xid(ci));
                    if length(tci) == 0
                        continue
                    end
                    if rlabel(ri) == rlabel(ci)
                        r_positive(tri,tci) = r_positive(tri,tci) + 1;
                    else
                        r_negtive(tri,tci) = r_negtive(tri,tci) - 1;
                    end
                end
            end
        end
        
        class_originals = zeros(3,nGroup); class_normals = zeros(3,nGroup);
        for iGroup = 1:nGroup
            if size(tempXs{iGroup},2)<3
                mean_norm_p(iGroup) = 1; mean_norm_d(iGroup) = 1;
            else
                class_originals(:,iGroup) =  mean(tempXs{iGroup},2);
                [U S V] = svd(tempXs{iGroup} - repmat(class_originals(:,iGroup),1,size(tempXs{iGroup},2)));
                class_normals(:,iGroup) = U(:,3);
                mean_norm_p(iGroup) = abs((Params.current - class_originals(:,iGroup))' * class_normals(:,iGroup));
                tempx_d{iGroup} = tempXs{iGroup} - repmat(Params.current,1,size(tempXs{iGroup},2));
                norm_tempx_d{iGroup} = sqrt(sum(tempx_d{iGroup}.*tempx_d{iGroup}));
                [value_d id_d] = sort(norm_tempx_d{iGroup});
                mean_norm_d(iGroup) = mean(value_d(1:min(6,length(value_d))));
            end
        end
        
        mean_norm_p = mean_norm_p / sum(mean_norm_p);
        mean_norm_d = mean_norm_d / sum(mean_norm_d);
        mean_norm = mean_norm_p + mean_norm_d;
        
        [value_g id_g] = min(mean_norm);
        norm_vector = class_normals(:,id_g)';
        break
    end
    nGroup = nGroup + 1;
end













