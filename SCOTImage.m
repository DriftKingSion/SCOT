%% Getting matrices from .mj2 or mp4. file

clear variables
Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec03.mj2');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec01.mj2');
Normal_frames=double(Normal_frames(:,:,:));
Dark_frames=double(Dark_frames(:,:,:));
% Dark_frames=imresize(Dark_frames,[480,852]);
% Normal_frames=squeeze(Normal_frames);

% Sizes of matrices
L_dark=size(Dark_frames,3);
L_normal=size(Normal_frames,3);
Height=size(Normal_frames,1);
Width=size(Normal_frames,2);

%% Getting source positions

N_sources=1;
Sources=zeros(N_sources,2);
minSize=50;

for i=1:N_sources
    Source=imbinarize(Normal_frames(:,:,i)/4095);    
    Source=bwareaopen(Source,minSize);
    s = regionprops(Source,'centroid');
    Sourcesmultiples = cat(1,s.Centroid);
    Sources(i,:) = Sourcesmultiples(1,:);
end

%% Creating an image with only the pixels in ROI

% Counts in electrons
CF=2.51;
Normal_frames_e=Normal_frames*CF;
Dark_frames_e=Dark_frames*CF;

% Substracting by mean of dark frames adn taking the mean
Dark_mean=mean(Dark_frames_e,3);
Normal_frame=mean(Normal_frames_e,3);
Darksubstracted_frame=Normal_frame-Dark_mean;

% One computation of distance for each source
% We attirbute each distance value to an indice which is stored in Distances_sources
for i=1:Height
    for j=1:Width
        D=pdist2([Sources(1)*0.00943,Sources(2)*0.012],[(j)*0.00943,(i)*0.012]); % Distances in mm
        if D<1 || D>2
            Normal_frame(i,j)=0;
            Darksubstracted_frame(i,j)=0;
        end
    end
end

%% Displaying the image

image(Darksubstracted_frame,'CDataMapping','scaled')
title('Display of average intensity')
colormap(hot)
colorbar