function plot_claasification_points(allpoints, knn, label)

id1 = find(label == 1);id2 = find(label == 2);id3 = find(label == 3);id4 = find(label == 4);id5 = find(label == 5);
knn1 = knn(id1); knn2 = knn(id2);knn3 = knn(id3);knn4 = knn(id4);knn5 = knn(id5);
h1= scatter3(allpoints(knn1,1),allpoints(knn1,2),allpoints(knn1,3),200,'.','MarkerEdgeColor',GS.CLASS_COLOR1);  hold on;
h2= scatter3(allpoints(knn2,1),allpoints(knn2,2),allpoints(knn2,3),200,'.','MarkerEdgeColor',GS.CLASS_COLOR2);  hold on;
h3= scatter3(allpoints(knn3,1),allpoints(knn3,2),allpoints(knn3,3),200,'.','MarkerEdgeColor',GS.CLASS_COLOR5); hold on;
h4= scatter3(allpoints(knn4,1),allpoints(knn4,2),allpoints(knn4,3),200,'.','MarkerEdgeColor',GS.CLASS_COLOR4); hold on;
h5= scatter3(allpoints(knn5,1),allpoints(knn5,2),allpoints(knn5,3),200,'.','MarkerEdgeColor',GS.CLASS_COLOR3); hold on;
delete(h1);delete(h2);delete(h3);delete(h4);delete(h5);

for iClass = 1:max(label)
    idi = find(label == iClass);knni = knn(idi);
    h1= scatter3(allpoints(knni,1),allpoints(knni,2),allpoints(knni,3),200,'.','MarkerEdgeColor',GS.PC_COLOR);  hold on;
    delete(h1);
end