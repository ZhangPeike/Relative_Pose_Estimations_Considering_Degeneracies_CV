function [ Image_Feature_Match_Error ] = Produce_Noise_for_Outlier( Image_Feature_Match,A )
%Make a right feature matches into error, turn inliers to be outliers
%Mechnism: For a 4xN input match, a set of same size random uniform number
%is produced, if a number is no less than 0.5, then the corresponding
%number is added with A, otherwise, -A.
    [M,N]=size(Image_Feature_Match);
    Image_Feature_Match_Error=zeros(M,N);
    if(M~=4)
        disp('Number of Match Row Error.');
    end
    N=size(Image_Feature_Match,2);
    Rand_Set=rand(M,N);
    for i=1:4
        for j=1:N
            if Rand_Set(i,j)>=0.5
                Image_Feature_Match_Error(i,j)=Image_Feature_Match(i,j)-A;
            else 
                Image_Feature_Match_Error(i,j)=Image_Feature_Match(i,j)+A;
            end
        end
    end
end

