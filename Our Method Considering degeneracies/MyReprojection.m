function [ ReprojectionPointHomo, ReprojectionPointInHomo ] = MyReprojection( P1,P2,X )
%将重建的三维空间点重投影到二维平面
%作者：张培科
%Input:相机矩阵，Camera Matrix 3*4，P1与P2，重建的三维点齐次坐标，n个点，矩阵4*n
%Output:二维点坐标，包括齐次的和非齐次的，x{1}齐次，x{2}非齐次，x{1}矩阵6*n，x{2}矩阵4*n
%可能存在比例，原始的二维图像点坐标以像素为尺度，而SFMedu获得的三维重建点的坐标尺度自己未知
%2016年3月
ReprojectionPointHomo(1:3,:)=P1*X;
ReprojectionPointHomo(4:6,:)=P2*X;
ReprojectionPointInHomo(1:2,:)=ReprojectionPointHomo(1:2,:)./ReprojectionPointHomo([3 3],:);
ReprojectionPointInHomo(3:4,:)=ReprojectionPointHomo(4:5,:)./ReprojectionPointHomo([6 6],:);
end

