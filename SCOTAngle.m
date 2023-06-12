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
    Sources(i,:) = cat(1,s.Centroid);
end

%% Computation of distances and K

% Counts in electrons
CF=2.51;
Normal_frames_e=Normal_frames*CF;
Dark_frames_e=Dark_frames*CF;

% Substracting by mean of dark frames
Dark_mean=mean(Dark_frames_e,3);
Darksubstracted_frames=Normal_frames_e-Dark_mean;

% Dividing all the possible angles by intervals
Nmb_of_points=60;
Ecart=2*pi/Nmb_of_points;
Angles_plot=linspace(-pi+Ecart/2,pi-Ecart/2,Nmb_of_points);

% Storing the values for every acquisition
Intn=zeros(Nmb_of_points,1);
Intd=zeros(Nmb_of_points,1);

% One computation of distance/angle for each source
% We attirbute each distance value to an indice which is stored in Distances_sources
Distances_sources=zeros(N_sources,Height-6,Width-6);
Angles_sources=zeros(N_sources,Height-6,Width-6);
for l=1:N_sources
    for i=1:Height-6
        for j=1:Width-6
            D=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+3)*0.00943,(i+3)*0.012]); % Distances in mm
            if D<1.4 || D>1.6
                Distances_sources(l,i,j)=0;
            else
                Distances_sources(l,i,j)=1;
            end
            A=angle(complex((j+3)-Sources(l,1),((i+3)-Sources(l,2))));
            for c=1:Nmb_of_points-1 
                if A>=Angles_plot(c) && A<Angles_plot(c+1)
                    Angles_sources(l,i,j)=c;
                end
            end
            if Angles_sources(l,i,j)==0
                Angles_sources(l,i,j)=Nmb_of_points;
            end
        end
    end
end

% Computation of variances and K
Compteur=zeros(Nmb_of_points,1);
for l=1:L_normal
    Normal_frame=Normal_frames_e(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Angle_source=Angles_sources(n,:,:);
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-6
        for j=1:Width-6
            if Distance_source(1,i,j)==1
                c=Angle_source(1,i,j);

                Mn=mean(Normal_frame(i:i+6,i:i+6),"all");
                Md=mean(Darksubstracted_frame(i:i+6,i:i+6),"all");

                Intn(c)=Intn(c)+Mn;
                Intd(c)=Intd(c)+Md;

                Compteur(c)=Compteur(c)+1;
            end
        end
    end
end

for c=1:Nmb_of_points
    Intn(c)=Intn(c)/Compteur(c);
    Intd(c)=Intd(c)/Compteur(c);
end

%% Plots

Intn2=[Intn ;Intn(1)];
Intd2=[Intd ;Intd(1)];
Angles_plot_2=[Angles_plot Angles_plot(1)];

polarplot(Angles_plot_2,Intd2)
title('Polar plot of the intensity averaged for all pixels between 1.4 and 1.6 mm from the source')