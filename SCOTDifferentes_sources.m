%% Getting matrices from .mj2 or mp4. file

clear variables
Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230523b\rec06.mj2');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230523b\rec02.mj2');
Normal_frames=double(Normal_frames(:,:,2:91));
Dark_frames=double(Dark_frames(:,:,1:50));
% Dark_frames=imresize(Dark_frames,[480,852]);
% Normal_frames=squeeze(Normal_frames);

% Sizes of matrices
L_dark=size(Dark_frames,3);
L_normal=size(Normal_frames,3);
Height=size(Normal_frames,1);
Width=size(Normal_frames,2);

%% Getting source positions

N_sources=6;
Sources=zeros(N_sources,2);
minSize=50;

for k=1:N_sources
    Source=imbinarize(Normal_frames(:,:,k)/4095);    
    Source=bwareaopen(Source,minSize);
    s = regionprops(Source,'centroid');
    Sources(k,:) = cat(1,s.Centroid);
end

%% Computation of distances and K

% Counts in electrons
CF=2.51;
Normal_frames=Normal_frames*CF;
Dark_frames=Dark_frames*CF;

% Substracting by mean of dark frames
Dark_mean=mean(Dark_frames,3);
Darksubstracted_frames=Normal_frames-Dark_mean;

% Dividing all the possible distances by intervals of size 0.05 mm
Nmb_of_points=51;
Dist_min=0.5;
Dist_max=3;
Ecart=(Dist_max-Dist_min)/(Nmb_of_points-1);
Distances_plot=linspace(Dist_min-Ecart,Dist_max+Ecart,Nmb_of_points+2);

% Storing the values for every acquisition
V_normal=zeros(Nmb_of_points,N_sources);
Intn=zeros(Nmb_of_points,N_sources);
Intd=zeros(Nmb_of_points,N_sources);
V_dark=zeros(Nmb_of_points,N_sources);
V_shot=zeros(Nmb_of_points,N_sources);
V_corr=zeros(Nmb_of_points,N_sources);
Speckle=zeros(Nmb_of_points,N_sources);
Speckle_corr=zeros(Nmb_of_points,N_sources);

% Size of the sliding pixel window
W_size=7;
RWindow=(W_size-1)/2; %Distance in pixels from center to the edge of the window

% Computing dark variance 
Mdv=0;
Compteur=0;
Mean_dark_V=zeros(Nmb_of_points,N_sources);
for l=1:L_dark
    Dark_frame=Dark_frames(:,:,l);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            Mdv=Mdv+var(Dark_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
            Compteur=Compteur+1;
        end
    end
end
Mdv=Mdv/Compteur;
Mean_dark_V(:)=Mdv;
clear('Dark_frames')

% One computation of distance for each source
% We attirbute each distance value to an indice which is stored in Distances_sources
Distances_sources=zeros(N_sources,Height-(W_size-1),Width-(W_size-1));
for l=1:N_sources
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            D=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+RWindow)*0.00943,(i+RWindow)*0.012]); % Distances in mm
            if D<Distances_plot(1)+Ecart/2 || D>Distances_plot(end)-Ecart/2
                Distances_sources(l,i,j)=nan;
            else
                c=1;
                while D>=Distances_plot(c)-Ecart/2
                    c=c+1;
                end
                Distances_sources(l,i,j)=c-2;
            end
        end
    end
end

% Computation of variances and K
Compteur=zeros(Nmb_of_points,N_sources);
for l=1:L_normal
    Normal_frame=Normal_frames(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            c=Distance_source(1,i,j);
            if ~isnan(c)
                Mn=mean(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Intn(c,n)=Intn(c,n)+Mn;
                Intd(c,n)=Intd(c,n)+Md;
                V_normal(c,n)=V_normal(c,n)+Vn;
                V_dark(c,n)=V_dark(c,n)+Vd;
                V_shot(c,n)=V_shot(c,n)+Md;
                V_corr(c,n)=V_corr(c,n)+(Vd-Mdv-Md);

                Compteur(c,n)=Compteur(c,n)+1;
            end
        end
    end
end

for n=1:N_sources
    for c=1:Nmb_of_points
        Intn(c,n)=Intn(c,n)/Compteur(c,n);
        Intd(c,n)=Intd(c,n)/Compteur(c,n);
        V_normal(c,n)=V_normal(c,n)/Compteur(c,n);
        V_dark(c,n)=V_dark(c,n)/Compteur(c,n);
        V_shot(c,n)=V_shot(c,n)/Compteur(c,n);
        V_corr(c,n)=V_corr(c,n)/Compteur(c,n);

        Speckle(c,n)=V_normal(c,n)/Intn(c,n)^2;
        Speckle_corr(c,n)=V_corr(c,n)/Intd(c,n)^2;
    end
end

%% Plots

t=tiledlayout(1,3);
title(t,'Computations for speckle contrast')

nexttile
plot(Distances_plot(2:end-1),log(Intn(:,1)),"+",Distances_plot(2:end-1),log(Intn(:,2)),"*",Distances_plot(2:end-1),log(Intn(:,3)),"x",Distances_plot(2:end-1),log(Intn(:,4)),"o",Distances_plot(2:end-1),log(Intn(:,5)),"x",Distances_plot(2:end-1),log(Intn(:,6)),"d")
xlabel('S-D separation (mm)')
ylabel('Log(Intensity) (e-)')
legend('Source1','Source2','Source3','Source4','Source5','Source6')
nexttile
plot(Distances_plot(2:end-1),Speckle(:,1),"+",Distances_plot(2:end-1),Speckle(:,2),"*",Distances_plot(2:end-1),Speckle(:,3),"x",Distances_plot(2:end-1),Speckle(:,4),"o",Distances_plot(2:end-1),Speckle(:,5),"x",Distances_plot(2:end-1),Speckle(:,6),"d")
xlabel('S-D separation (mm)')
ylabel('K^{2}')
legend('Source1','Source2','Source3','Source4','Source5','Source6')
nexttile
plot(Distances_plot(2:end-1),Speckle_corr(:,1),"+",Distances_plot(2:end-1),Speckle_corr(:,2),"*",Distances_plot(2:end-1),Speckle_corr(:,3),"x",Distances_plot(2:end-1),Speckle_corr(:,4),"o",Distances_plot(2:end-1),Speckle_corr(:,5),"x",Distances_plot(2:end-1),Speckle_corr(:,6),"d")
xlabel('S-D separation (mm)')
ylabel('K_{c}^{2}')
legend('Source1','Source2','Source3','Source4','Source5','Source6')