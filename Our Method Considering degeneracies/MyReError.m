function [ ErrorSum ] = MyReError( matches,reprojection )
%Input 输入：matches，原图像提取的特征点对坐标，n对点，矩阵4*n，reprojection为重投影生成的点对坐标，格式同matches
%Output输出：原始点与对应重投影点的欧氏距离的平方和
%不可作为优化指标
ErrorSum=(matches(1,1)-reprojection(1,1))^2+(matches(2,1)-reprojection(2,1))^2+(matches(3,1)-reprojection(3,1))^2+(matches(4,1)-reprojection(4,1))^2;
    for i=2:size(matches,2)
        ErrorSum=ErrorSum+(matches(1,i)-reprojection(1,i))^2+(matches(2,i)-reprojection(2,i))^2+(matches(3,i)-reprojection(3,i))^2+(matches(4,i)-reprojection(4,i))^2;
    end
end

