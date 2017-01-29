function [ Right_Rt, Qt, D, M] = Universial_Relative_PoseV3( K1,K2,matches,tF,tH)
% A Universial Solution to Relative Pose Problem Considering motion degeneracy and structrue degeneracy
% Author:Zhang Peike
% Date:2016.Sept.10
% Version: 3
% Input:tF, threshold for the fundamental matrix; tH, threshold for the homography matrix
% Date: 2017-01-16 18:49:31
% Detail: Using RANSAC to be robust.
%% Fundamental matrix
% Robust estimating F
%[ F,inliersF ] = FRANSACScore( matches,t,T,N );
%Reprojection error threshold 
ERROR=500;
%[F,inliersF] = ransacfitfundmatrix(matches(1:2,:),matches(3:4,:),tF);
%% Homography matrix
% Estimating H Robustly with linear method
%[ H,inliers ] = HRANSAC( matches,t,T,N )
%H = vgg_H_from_x_lin(matches(1:2,inliersF),matches(3:4,inliersF));
[H,inliersH] = ransacfithomography_Est_with_Normalizing_Err_original(matches(1:2,:),matches(3:4,:),tH);
%[H,inliersH] = ransacfithomography_vgg_linear(matches(1:2,:),matches(3:4,:),tH);
%disp('H is done!');
%disp(H);
% Transfering to General matrix and recover pose
%G = K2\H*K1;
Qt = true;
[ RtH,isPureRotation ] = HtoRtV3( H,K1,K2,matches(:,inliersH) );
    if  isPureRotation
        Right_Rt=RtH;
        Qt=true;
        D=0;
        M=3;
        return;
    else
        [F,inliersF,D] = ransacfitfundmatrixDistanceWithoutNormalization(matches(1:2,:),matches(3:4,:),tF);
        % Transfering to E,then Rt
        if D==0
            E=K2'*F*K1;
            %disp('F&E is done!');
            %disp(E);
            % No using SVD to get Pose from Essential
            %[ RightRt ] = EtoRtN( E,K1,K2,matches )
            RtE = EtoRtV3( E,K1,K2,matches(:,inliersF) );
            % Reprojecting and computing error
            %[ Optmised_X ] = TriangulationOpt( K1,K2,R1,R2,t1,t2,matches,imsize)
            %X_Triangualed_E = TriangulationOptV2(K1,K2,eye(3),RtE(:,1:3),zeros(3,1),RtE(:,4),matches(:,inliersF));
            P1=K1*[eye(3) zeros(3,1)];
            P2=K2*RtE;
            X_Triangualed_E = Triangulation(P1,P2,matches(:,inliersF));
            %X_Triangualed_E_Inhomo=X_Triangualed_E(1:3,:)./X_Triangualed_E([4,4,4],:);
            [ ~, ReprojectionPointInHomoE ] = MyReprojection( P1,P2,X_Triangualed_E );
            %{
            figure;
            plot3(X_Triangualed_F_Inhomo(1,:),X_Triangualed_F_Inhomo(2,:),X_Triangualed_F_Inhomo(3,:),'.g');
            grid on
            title('F');
            %}
            [Error_List_E,Error_Sum_E]= MyReErrorIndi( matches(:,inliersF),ReprojectionPointInHomoE );   
            %X_Triangualed_H = Triangulation(K1*[eye(3) zeros(3,1)],K2*RtH,matches(:,inliersF));
            %X_Triangualed_H = TriangulationOptV2(K1,K2,eye(3),RtH(:,1:3),zeros(3,1),RtH(:,4),matches(:,inliersH));
            P2=K2*RtH;
            X_Triangualed_H = Triangulation(P1,P2,matches(:,inliersH));
            %X_Triangualed_H_Inhomo = X_Triangualed_H(1:3,:)./X_Triangualed_H([4,4,4],:);
            [ ~, ReprojectionPointInHomoH ] = MyReprojection( P1,P2,X_Triangualed_H );
            %{
            figure;
            plot3(X_Triangualed_H_Inhomo(1,:),X_Triangualed_H_Inhomo(2,:),X_Triangualed_H_Inhomo(3,:),'.g');
            grid on
            title('H');
            %}
            [Error_List_H,Error_Sum_H] = MyReErrorIndi( matches(:,inliersH),ReprojectionPointInHomoH );
            meanErrorE=Error_Sum_E/length(inliersF);
            meanErrorH=Error_Sum_H/length(inliersH);        
            medianErrorE=median(Error_List_E);
            medianErrorH=median(Error_List_H);        
            Qt=true;
            if(medianErrorE>ERROR&&medianErrorH>ERROR)
                Qt=false;
            end
            % 2017-01-21 20:09:43
            % weight for the model function of number of inliers and the
            % reprojection error
            w_H_1=0;
            w_H_2=0;
            w_E_1=0;
            w_E_2=0;
            N_H=length(inliersH);
            N_E=length(inliersF);
            if N_H<0.25*N_E
                w_H_1=0.25;
            else if N_H<0.75*N_E && N_H>=0.25*N_E
                    w_H_1=0.5;
                else if N_H>=0.75*N_E && N_H<=1.25*N_E
                        w_H_1=1;
                    else if N_H>1.25*N_E && N_H <=2*N_E
                            w_H_1=2;
                        else if N_H>2*N_E
                                w_H_1=4;
                            end
                        end
                    end
                end
            end
            if N_E<0.25*N_H
                w_E_1=0.25;
            else if N_E>=0.25*N_H && N_E<0.75*N_H
                    w_E_1=0.5;
                else if N_E>=0.75*N_H && N_E<=1.25*N_H
                        w_E_1=1;
                    else if N_E>1.25*N_H && N_E<=2*N_H
                            w_E_1=2;
                        else if N_E>2*N_H
                                w_E_1=4;
                            end
                        end
                    end
                end
            end
            Max_Err=max(medianErrorE,medianErrorH);
            w_H_2=exp(1-(medianErrorH/Max_Err));
            w_E_2=exp(1-(medianErrorE/Max_Err));
            w_H=w_H_1*w_H_2;
            w_E=w_E_1*w_E_2;
            if w_H<w_E
                Right_Rt=RtE;
                disp('E is chosen.');
                M=1;
            else if w_E<w_H
                    Right_Rt=RtH;
                    disp('H is chosen.');
                    M=2;
                else
                    M=0;
                    disp('Scores of H and E are equal');
                end
            end
        else
            Right_Rt=zeros(3,4);
            Qt=false;
            M=0;
        end
    end
end
