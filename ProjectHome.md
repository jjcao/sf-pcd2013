# Point cloud normal estimation via low-rank subspace clustering #
Jie Zhang, Junjie Cao, Xiuping Liu, Jun Wang, Jian Liu, Xiquan Shi<br />
Computer & Graphics (SMI 2013)

<img src='http://sf-pcd.googlecode.com/svn/trunk/write_up/figures/teaser.png' width='800' title='Teaser' /><br />
Reconstructed normals of two plans with shallow angle by <Li et al. 10> (left), <Boulch et al. 12> (middle) and our algorithm (right), respectively. The points are sampled uniformly in the top row and non-uniformly in the bottom row.<br />

## Motivation ##
  * Estimating surface normals from a noisy point cloud benefits many applications in computer graphics, geometry processing and reverse engineering.
  * In point cloud models noise, sharp features, outliers and non-uniform distribution are ubiquitous. These problems make it difficult to estimate normals accurately. When near sharp features even the best normal estimation algorithm still lead to some errors. While subspace clustering provide an effective method for subneighborhood selection which can assist in normal estimation.

## Abstract ##
In this paper, we present a robust normal estimation algorithm based on the low-rank subspace clustering technique. The main idea is based on the observation that compared with the points around sharp features, it is relatively easier to obtain accurate normals for the points within smooth regions. The points around sharp features and smooth regions are identified by covariance analysis of their neighborhoods. The neighborhood of a point in a smooth region can be well approximated by a plane. For a point around sharp features, some of its neighbors may be in smooth regions. These neighbor points' normals are estimated by principal component analysis, and used as prior knowledge to carry out neighborhood clustering. An unsupervised learning process is designed to represent the prior knowledge as a guiding matrix. Then we segment the anisotropic neighborhood into several isotropic neighborhoods by low-rank subspace clustering with the guiding matrix, and identify a consistent subneighborhood for the current point. Hence the normal of the current point near sharp features is estimated as the normal of a plane fitting the consistent subneighborhood. Our method is capable of estimating normals accurately even in the presence of noise and anisotropic samplings, while preserving sharp features within the original point data. We demonstrate the effectiveness and robustness of the proposed method on a variety of examples.

## Keyword ##
normal estimation; sharp feature preserving; low-rank representation; subspace clustering

## Results ##
<img src='http://sf-pcd.googlecode.com/svn/trunk/write_up/figures/fandisk_mark.jpg' width='800' title='fandisk' /><br />
Comparison for Fandisk model with 50% noise. From top to bottom row are the results rendering using surfels, the visualization of bad points and top view of the generated normals, respectively. Left to right are the results of <Hoppe et al. 92>,<Li et al. 10>, <Boulch et al. 12> and our algorithm, respectively.<br />

<img src='http://sf-pcd.googlecode.com/svn/trunk/write_up/figures/shutter.png' width='800' title='shutter' /><br />
Normal estimation for the Shutter model. Left to right are the input raw scan model with 291.2k points, the results of <Hoppe et al. 92> and our algorithm, respectively.<br />

## Applications ##
### Sharp Feature Detection ###
<img src='http://sf-pcd.googlecode.com/svn/trunk/write_up/figures/feature.png' width='400' title='feature' /><br />
The sharp feature detect results by <Park et al. 12> (top) and our method (bottom). Left: Smooth-feature model with 10% Gaussian noise. Right: Smooth-feature model with 20% Gaussian noise.<br />

## Bibtex ##
@article{DBLP:journals/cg/ZhangCLWLS13, Author = {Jie Zhang and unjie Cao and Xiuping Liu and Jun Wang and Jian Liu and Xiquan Shi}, Title = {Point cloud normal estimation via low-rank subspace clustering}, Journal = {Computers {\&} Graphics}, volume  = {37}, number = {6}, year = {2013}, pages = {697-706}}<br />
Jie Zhang, Junjie Cao, Xiuping Liu, Jun Wang, Jian Liu, Xiquan Shi. Point cloud normal estimation via low-rank subspace clustering. Computer & Graphics (SMI 2013), 2013, 37(6): 697-706.<br />

## Main References ##
  1. <Hoppe et al. 92> H. Hoppe, T. DeRose, T. Duchamp, J. A. McDonald,W. Stuetzle, Surface reconstruction from unorganized points, in: SIGGRAPH, 1992, pp. 71–78.
  1. <Li et al. 10> B. Li, R. Schnabel, R. Klein, Z.-Q. Cheng, G. Dang, S. Jin, Robust normal estimation for point clouds with sharp features, Computers & Graphics34 (2) (2010) 94–106.
  1. <Boulch et al. 12> A. Boulch, R. Marlet, Fast and robust normal estimation for point clouds with sharp features, Comput. Graph. Forum  31 (5) (2012) 1765–1774.
  1. <Park et al. 12> M. K. Park, S. J. Lee, K. H. Lee, Multi-scale tensor voting for feature extraction from unstructured point clouds, Graphical Models 74 (4) (2012) 197–208.