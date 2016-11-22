function [ Right_Rt ] = Universial_Relative_PoseV2s( K1,K2,matches,threshold )
% Paper: A Universial Solution to Relative Pose Problem
% Considering motion degeneracy and structrue degeneracy
% Author:Zhang Peike
% Date:2016.Sept.10
% Detail: Inliers are determined by the F matrix estimating
% Remaining task: Estimating Homography nonlinearly 
% imsize: Width Height, opposite to the size function!!!
%% Fundamental matrix
% Robust estimating F
%[ F,inliersF ] = FRANSACScore( matches,t,T,N )
[F,inliersF]=ransacfitfundmatrix(matches(1:2,:),matches(3:4,:),threshold);
% Transfering to E,then Rt
E=K2'*F*K1;
%disp('F&E is done!');
%disp(E);
%% Homography matrix
% Estimating H linearly
%[ H,inliers ] = HRANSAC( matches,t,T,N )
H = vgg_H_from_x_lin(matches(1:2,inliersF),matches(3:4,inliersF));
%disp('H is done!');
%disp(H);
% Transfering to General matrix and recover pose
%G = K2\H*K1;
[ RtH,isPureRotation ] = Rt_HOptV2( H,K1,K2,matches(:,inliersF) );
if  isPureRotation
    Right_Rt=RtH;
else
    % No using SVD to get Pose from Essential
    %[ RightRt ] = EtoRtN( E,K1,K2,matches )
    RtF = EtoRtV3( E,K1,K2,matches(:,inliersF) );
    % Reprojecting and computing error
    %[ Optmised_X ] = TriangulationOpt( K1,K2,R1,R2,t1,t2,matches,imsize)
    X_Triangualed_F = TriangulationOptV2(K1,K2,eye(3),RtF(:,1:3),zeros(3,1),RtF(:,4),matches(:,inliersF));
    X_Triangualed_F_Inhomo = X_Triangualed_F(1:3,:)./X_Triangualed_F([4,4,4],:);
    [ ~, ReprojectionPointInHomoF ] = MyReprojection( K1*[eye(3) zeros(3,1)],K2*RtF,X_Triangualed_F );
    %{
    figure;
    plot3(X_Triangualed_F_Inhomo(1,:),X_Triangualed_F_Inhomo(2,:),X_Triangualed_F_Inhomo(3,:),'.g');
    grid on
    title('F');
    %}
    ErrorSum_F = MyReError( matches(:,inliersF),ReprojectionPointInHomoF );   
    %X_Triangualed_H = Triangulation(K1*[eye(3) zeros(3,1)],K2*RtH,matches(:,inliersF));
    X_Triangualed_H = TriangulationOptV2(K1,K2,eye(3),RtH(:,1:3),zeros(3,1),RtH(:,4),matches(:,inliersF));
    X_Triangualed_H_Inhomo = X_Triangualed_H(1:3,:)./X_Triangualed_H([4,4,4],:);
    [ ~, ReprojectionPointInHomoH ] = MyReprojection( K1*[eye(3) zeros(3,1)],K2*RtH,X_Triangualed_H );
    %{
    figure;
    plot3(X_Triangualed_H_Inhomo(1,:),X_Triangualed_H_Inhomo(2,:),X_Triangualed_H_Inhomo(3,:),'.g');
    grid on
    title('H');
    %}
    ErrorSum_H = MyReError( matches(:,inliersF),ReprojectionPointInHomoH );
    if ErrorSum_H<ErrorSum_F
        Right_Rt=RtH;
        disp('Chosen H');
    else
        Right_Rt=RtF;
        disp('Chosen F');
    end  
end
end

