%% Getting matrices from .mj2 or mp4. file

Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec03.mj2');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec01.mj2');
Normal_frames=double(Normal_frames(:,:,1:500));
Dark_frames=double(Dark_frames(:,:,1:200));

% Sizes of matrices
L_dark=size(Dark_frames,3);
L_normal=size(Normal_frames,3);
Height=size(Normal_frames,1);
Width=size(Normal_frames,2);

%% Getting source positions

N_sources=1;
Sources=zeros(N_sources,2);
minSize = 50;

for k=1:N_sources
    Source=imbinarize(Normal_frames(:,:,k)/4095);    
    Source=bwareaopen(Source,minSize);
    s = regionprops(Source,'centroid');
    Sources(k,:) = cat(1,s.Centroid);
end

%% Computation of distances

% One computation for each source
Distances_sources_2=zeros(N_sources,Height-6,Width-6);
for l=1:N_sources
    for i=1:Height-6
        for j=1:Width-6
            Distances_sources_2(l,i,j)=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+3)*0.00943,(i+3)*0.012]);  % Distance in mm
        end
    end
end

% Creation of the complete matrix for every frame
Distances_sources=Distances_sources_2;
for l=1:L_normal/N_sources-1
    Distances_sources=cat(1,Distances_sources,Distances_sources_2);
end

%% Computation of K

% Counts in electrons
CF=2.51;
Normal_frames_e=Normal_frames*CF;
Dark_frames_e=Dark_frames*CF;

% Substracting by mean of dark frames
Dark_mean=mean(Dark_frames_e,3);
Darksubstracted_frames=Normal_frames_e-Dark_mean;

% Computing dark variance
Mean_dark_V=zeros(L_normal,(Height-6),(Width-6)); 
Var_dark=zeros(L_normal,(Height-6),(Width-6));
for l=1:L_dark
    for i=1:Height-6
        for j=1:Width-6
        Dark_frame=Dark_frames(:,:,l);
        Var_dark(l,i,j)=var(Dark_frame(i:i+6,j:j+6),0,"all");
        end
    end
end
Mean_dark_V(:,:,:)=mean(Var_dark,"all");

% Computing other variances

% Storing the values for every acquisition
V_normal=zeros(L_normal,Height-6,Width-6);
Int=zeros(L_normal,Height-6,Width-6);
V_dark=zeros(L_normal,Height-6,Width-6);
V_shot=zeros(L_normal,Height-6,Width-6);
V_corr=zeros(L_normal,Height-6,Width-6);
Speckle=zeros(L_normal,Height-6,Width-6);
Speckle_corr=zeros(L_normal,Height-6,Width-6);

for l=1:L_normal
    Normal_frame=Normal_frames_e(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    for i=1:Height-6
        for j=1:Width-6
            if Distances_sources(l,i,j)>1 && Distances_sources(l,i,j)<1.5
            Int(l,i,j)=mean(Darksubstracted_frame(i:i+6,j:j+6),"all");
            V_normal(l,i,j)=var(Normal_frame(i:i+6,j:j+6),0,"all");
            V_dark(l,i,j)=var(Darksubstracted_frame(i:i+6,j:j+6),0,"all");
            V_shot(l,i,j)=Int(l,i,j);
            V_corr(l,i,j)=V_dark(l,i,j)-Mean_dark_V(l,i,j)-V_shot(l,i,j);

            Speckle(l,i,j)=sqrt(V_normal(l,i,j))/(Int(l,i,j)^2);
            Speckle_corr(l,i,j)=sqrt((V_corr(l,i,j))/(Int(l,i,j)^2));
            end
            if Distances_sources(l,i,j)<=1 || Distances_sources(l,i,j)>=1.5
            V_normal(l,i,j)=nan;
            Int(l,i,j)=nan;
            V_dark(l,i,j)=nan;
            V_shot(l,i,j)=nan;
            V_corr(l,i,j)=nan;
            Speckle(l,i,j)=nan;
            Speckle_corr(l,i,j)=nan;
            end
        end
    end
end 

%% Putting the data in 2D and unit16

Distances_sources=reshape(Distances_sources,L_normal*(Height-6)*(Width-6),1);
Int=reshape(Int,L_normal*(Height-6)*(Width-6),1);
V_dark=reshape(V_dark,L_normal*(Height-6)*(Width-6),1);
V_shot=reshape(V_shot,L_normal*(Height-6)*(Width-6),1);
V_corr=reshape(V_corr,L_normal*(Height-6)*(Width-6),1);
Speckle=reshape(Speckle,L_normal*(Height-6)*(Width-6),1);
Speckle_corr=reshape(Speckle_corr,L_normal*(Height-6)*(Width-6),1);
Mean_dark_V=reshape(Mean_dark_V,L_normal*(Height-6)*(Width-6),1);

%% Plots

t=tiledlayout(2,2);
title(t,'Computations for speckle contrast')
nexttile
plot(Distances_sources,Int,"+")
xlabel('S-D separation (mm)')
ylabel('Intensity')
nexttile
plot(Distances_sources,V_dark,"*",Distances_sources,V_shot,"x",Distances_sources,V_corr,"o")
xlabel('S-D separation (mm)')
ylabel('Variance')
legend("Variance of dark substracted intensity","Shot noise variance","Corrected variance")
nexttile
plot(Int,Mean_dark_V,Int,V_dark,"*",Int,V_shot,"x",Int,V_corr,"o")
xlabel('Intensity')
ylabel('Variance')
legend("Dark Variance","Variance of dark substracted intensity","Shot noise variance","Corrected variance")
nexttile
plot(Distances_sources,Speckle,"+",Distances_sources,Speckle_corr,"*")
xlabel('S-D separation (mm)')
ylabel('K')
legend("Uncorrected speckle contrast","Corrected speckle contrast")
