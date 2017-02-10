Relative_Pose_Estimations_Considering_Degeneracies_CV
This repo including functions and dataset for testing traditional and our novel method in estimating camera relative pose.
For the real image data, we have managed to obtain point feature correspondence by typical SIFT algorithm and get the camera calibrated. 
Results (matches and K matrix) can be used directly to test the relative pose estimating algorithms.
The model estimations are all run in RANSAC loop to avoid possible outliers.
The functions are written in Matlab (R2015 or later) for windows.
We have used RANSAC scheme and thirdparty software(Pix4D, Multiple View Environment) and we are adding our method into some noteable visual-SLAM to improve its tracking (pose esitmation) performance.
We have managed to reconstruct planar structure with the proposed pipeline embedded into the initialization step of the SFM part of MVE, please see the repo 3D-reconstruction-system for the soucre code and built exetuble.
Feel free to contact, zhangpeike@csu.edu.cn
Author: Peike Zhang,Yuanxin Wu, Qi Cai and Danping Zou
Some References
B. K. Horn, "Recovering baseline and orientation from essential matrix," J. Opt. Soc. Am, vol. 110, 1990.
R. Hartley and A. Zisserman, Multiple view geometry in computer vision: Cambridge university press, 2003.
D. Nister, "An efficient solution to the five-point relative pose problem," in Computer Vision and Pattern Recognition, 2003. Proceedings. 2003 IEEE Computer Society Conference on, 2003, pp. II-195-202 vol.2.
O. D. Faugeras and F. Lustman, "Motion and structure from motion in a piecewise planar environment," International Journal of Pattern Recognition and Artificial Intelligence, vol. 2, pp. 485-508, 1988.
G. H. Golub and C. F. Van Loan, Matrix Computations: Johns Hopkins University Press, 2012



