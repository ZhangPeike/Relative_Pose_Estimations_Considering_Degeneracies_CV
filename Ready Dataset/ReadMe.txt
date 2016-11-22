These Matlab mat format files consist of all kinds of data which is necessary to test the relative pose estimating algorithms.

The cells 'Mathch' or 'Frame' include several 4 * n arrays, with every column of the array meaning a feature point correspondence across two images (x1,y1,x2,y2)... .
 
The K1£¬K2 or K is the camera(s)' internal parameter matrix.

The synthetic data file contains four configurations: 1.General structure and general motion, 2. Planar structure and general motion, 3.General structure and pure translation and 4.Planar structure and pure translation.
The ground truth is also contained.The attitude is represented in quaternion and stored in the cell REAL_Q. The RealC_GG and RealC_DG data means the displacement of the camera center in configuration 1 and 2 respectively.

In the following files, feature point matches have been obtained by using SIFT 
algorithms. The camera internal matrix is acquired by the Matlab camera calibration toolbox.

In the file Doodle_Wall_Data.mat, note the image scene represents a typical planar structure.

In the file Planar_Structure_Sequence_Data.mat, note that the real movement of the camera corresponding to sequencing inter-frame is£ºgeneral motion, pure rotation, pure rotation and general motion.

In the file General_Structure_Sequence_Data.mat, note that the real movement of the camera corresponding to sequencing inter-frame is£ºpure rotation, pure rotation, pure rotation, pure rotation and general motion.

In the file Pure_Rotation_Sequence_Data.mat, note that the real movement of the camera corresponding to sequencing inter-frame is all pure rotation.

END






