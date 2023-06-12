%% Getting matrices from .mj2 or mp4. file

clear variables
Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec03.mj2');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230530\rec01.mj2');
Normal_frames=double(Normal_frames(:,:,1:1000));
Dark_frames=double(Dark_frames(:,:,1:200));
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

for k=1:N_sources
    Source=imbinarize(Normal_frames(:,:,k)/4095);    
    Source=bwareaopen(Source,minSize);
    s = regionprops(Source,'centroid');
    Sources(k,:) = cat(1,s.Centroid);
end

%% Computation of distances and K

% Counts in electrons
CF=2.51;
Normal_frames_e=Normal_frames*CF;
Dark_frames_e=Dark_frames*CF;

% Substracting by mean of dark frames
Dark_mean=mean(Dark_frames_e,3);
Darksubstracted_frames=Normal_frames_e-Dark_mean;

% Dividing all the possible distances by intervals of size 0.05 mm
Nmb_of_points=51;
Dist_min=0.5;
Dist_max=3;
Ecart_d=(Dist_max-Dist_min)/(Nmb_of_points);
Distances_plot=linspace(Dist_min+Ecart_d/2,Dist_max-Ecart_d/2,Nmb_of_points);

% Storing the values for every acquisition and for every number of frames
V_normal=zeros(Nmb_of_points,3);
Intn=zeros(Nmb_of_points,3);
Intd=zeros(Nmb_of_points,3);
V_dark=zeros(Nmb_of_points,3);
V_shot=zeros(Nmb_of_points,3);
V_corr=zeros(Nmb_of_points,3);
Speckle=zeros(Nmb_of_points,3);
Speckle_corr=zeros(Nmb_of_points,3);

% Size of the sliding pixel window
W_size=7;
RWindow=(W_size-1)/2; %Distance in pixels from center to the edge of the window

% Computing dark variance 
Mdv=0;
Compteur=0;
Mean_dark_V=zeros(Nmb_of_points,1);
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

% One computation of distance for each source
% We attirbute each distance value to an indice which is stored in Distances_sources
Distances_sources=zeros(N_sources,Height-6,Width-6);
for l=1:N_sources
    for i=1:Height-6
        for j=1:Width-6
            D=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+3)*0.00943,(i+3)*0.012]); % Distances in mm
            if D<Distances_plot(1)-(Ecart_d/2) || D>Distances_plot(end)+(Ecart_d/2)
                Distances_sources(l,i,j)=nan;
            else
                for c=1:Nmb_of_points 
                    if D>=Distances_plot(c)-(Ecart_d/2) && D<Distances_plot(c)+(Ecart_d/2)
                        Distances_sources(l,i,j)=c;
                    end
                end
            end
        end
    end
end

% Number of frames we want to test
N_frames=[5,50,500];
Compteur=zeros(Nmb_of_points,3);
Sigma=zeros(Nmb_of_points,3);

%% First group of frames
z=1;
N=N_frames(z);
Frame_number=randperm(L_normal,N);

% Computation of variances and K
for o=1:N
    l=Frame_number(o);
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
                Intn(c,z)=Intn(c,z)+Mn;
                Intd(c,z)=Intd(c,z)+Md;
                V_normal(c,z)=V_normal(c,z)+Vn;
                V_dark(c,z)=V_dark(c,z)+Vd;
                V_shot(c,z)=V_shot(c,z)+Md;
                V_corr(c,z)=V_corr(c,z)+(Vd-Mdv-Md);

                Compteur(c,z)=Compteur(c,z)+1;
            end
        end
    end
end

Intn(:,z)=Intn(:,z)./Compteur(:,z);
Intd(:,z)=Intd(:,z)./Compteur(:,z);
V_normal(:,z)=V_normal(:,z)./Compteur(:,z);
V_dark(:,z)=V_dark(:,z)./Compteur(:,z);
V_shot(:,z)=V_shot(:,z)./Compteur(:,z);
V_corr(:,z)=V_corr(:,z)./Compteur(:,z);

Speckle(:,z)=V_normal(:,z)./Intn(:,z).^2;
Speckle_corr(:,z)=V_corr(:,z)./Intd(:,z).^2;

% Computaion of the standard deviation of K for each value

for l=1:N
    Normal_frame=Normal_frames(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            c=Distance_source(1,i,j);
            if ~isnan(c)
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                ecart_type=(((Vd-Mdv-Md)/Md^2)-Speckle_corr(c,z))^2;
                Sigma(c,z)=Sigma(c,z)+ecart_type;
            end
        end
    end
end

%% Second group of frames
z=2;
N=N_frames(z);
Frame_number=randperm(L_normal,N);

% Computation of variances and K
for o=1:N
    l=Frame_number(o);
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
                Intn(c,z)=Intn(c,z)+Mn;
                Intd(c,z)=Intd(c,z)+Md;
                V_normal(c,z)=V_normal(c,z)+Vn;
                V_dark(c,z)=V_dark(c,z)+Vd;
                V_shot(c,z)=V_shot(c,z)+Md;
                V_corr(c,z)=V_corr(c,z)+(Vd-Mdv-Md);

                Compteur(c,z)=Compteur(c,z)+1;
            end
        end
    end
end

Intn(:,z)=Intn(:,z)./Compteur(:,z);
Intd(:,z)=Intd(:,z)./Compteur(:,z);
V_normal(:,z)=V_normal(:,z)./Compteur(:,z);
V_dark(:,z)=V_dark(:,z)./Compteur(:,z);
V_shot(:,z)=V_shot(:,z)./Compteur(:,z);
V_corr(:,z)=V_corr(:,z)./Compteur(:,z);

Speckle(:,z)=V_normal(:,z)./Intn(:,z).^2;
Speckle_corr(:,z)=V_corr(:,z)./Intd(:,z).^2;

% Computaion of the standard deviation of K for each value
for l=1:N
    Normal_frame=Normal_frames(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            c=Distance_source(1,i,j);
            if ~isnan(c)
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                ecart_type=(((Vd-Mdv-Md)/Md^2)-Speckle_corr(c,z))^2;
                Sigma(c,z)=Sigma(c,z)+ecart_type;
            end
        end
    end
end

%% Last group of frames

z=3;
N=N_frames(z);
Frame_number=randperm(L_normal,N);

% Computation of variances and K
for o=1:N
    l=Frame_number(o);
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
                Intn(c,z)=Intn(c,z)+Mn;
                Intd(c,z)=Intd(c,z)+Md;
                V_normal(c,z)=V_normal(c,z)+Vn;
                V_dark(c,z)=V_dark(c,z)+Vd;
                V_shot(c,z)=V_shot(c,z)+Md;
                V_corr(c,z)=V_corr(c,z)+(Vd-Mdv-Md);

                Compteur(c,z)=Compteur(c,z)+1;
            end
        end
    end
end

Intn(:,z)=Intn(:,z)./Compteur(:,z);
Intd(:,z)=Intd(:,z)./Compteur(:,z);
V_normal(:,z)=V_normal(:,z)./Compteur(:,z);
V_dark(:,z)=V_dark(:,z)./Compteur(:,z);
V_shot(:,z)=V_shot(:,z)./Compteur(:,z);
V_corr(:,z)=V_corr(:,z)./Compteur(:,z);

Speckle(:,z)=V_normal(:,z)./Intn(:,z).^2;
Speckle_corr(:,z)=V_corr(:,z)./Intd(:,z).^2;

% Computaion of the standard deviation of K for each valu
for l=1:N
    Normal_frame=Normal_frames(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            c=Distance_source(1,i,j);
            if ~isnan(c)
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                ecart_type=(((Vd-Mdv-Md)/Md^2)-Speckle_corr(c,z))^2;
                Sigma(c,z)=Sigma(c,z)+ecart_type;
            end
        end
    end
end

%% Compute the final matrix for standard deviation

Sigma_2=sqrt(Sigma./Compteur);
Sigma_plot=zeros(Nmb_of_points,3);
for l=1:Nmb_of_points
    n=mod(l,3)+1;
    Sigma_plot(l,n)=Sigma_2(l,n);
end

%% Plots
Distances_plot_2=[Distances_plot.' Distances_plot.' Distances_plot.'];
errorbar(Distances_plot_2,Speckle_corr,Sigma_plot,"o")
xlabel('S-D separation (mm)')
ylabel('K^2')
title('Corrected speckle contrast with errorbars')
legend([num2str(N_frames(1)) 'images'],[num2str(N_frames(2)) 'images'],[num2str(N_frames(3)) 'images'])