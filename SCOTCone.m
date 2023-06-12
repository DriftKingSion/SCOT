%% Getting matrices from .mj2 or mp4. file

clear variables
Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230607\rec02.mj2');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230607\rec01.mj2');
Normal_frames=double(Normal_frames(:,:,1:200));
Dark_frames=double(Dark_frames(:,:,1:50));
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
minSize=80;

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
Dist_min=1;
Dist_max=5;
Ecart_d=(Dist_max-Dist_min)/(Nmb_of_points);
Distances_plot=linspace(Dist_min+Ecart_d/2,Dist_max-Ecart_d/2,Nmb_of_points);

% Accepted angles for the computation (in rad)
Angle_min=-10*2*pi/360;
Angle_max=10*2*pi/360;

% Storing the values for every acquisition
V_normal=zeros(Nmb_of_points,1);
Intn=zeros(Nmb_of_points,1);
Intd=zeros(Nmb_of_points,1);
V_dark=zeros(Nmb_of_points,1);
V_shot=zeros(Nmb_of_points,1);
V_corr=zeros(Nmb_of_points,1);
Speckle=zeros(Nmb_of_points,1);
Speckle_corr=zeros(Nmb_of_points,1);

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
Angles_sources=zeros(N_sources,Height-6,Width-6);
for l=1:N_sources
    for i=1:Height-6
        for j=1:Width-6
            D=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+3)*0.00943,(i+3)*0.012]); % Distances in mm
            if D<Distances_plot(1)-(Ecart_d/2) || D>Distances_plot(end)+(Ecart_d/2)
                Distances_sources(l,i,j)=nan;
            else
                for d=1:Nmb_of_points 
                    if D>=Distances_plot(d)-(Ecart_d/2) && D<Distances_plot(d)+(Ecart_d/2)
                        Distances_sources(l,i,j)=d;
                    end
                end
            end
            A=angle(complex((j+3)-Sources(l,1),((i+3)-Sources(l,2))));
            if A<Angle_min || A>Angle_max
                Angles_sources(l,i,j)=nan;
            else
                Angles_sources(l,i,j)=1;
            end
        end
    end
end

% Computation of variances and K
Compteur=zeros(Nmb_of_points,1);
for l=1:L_normal
    Normal_frame=Normal_frames(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;
    Distance_source=Distances_sources(n,:,:);
    Angle_source=Angles_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            a=Angle_source(1,i,j);
            d=Distance_source(1,i,j);
            if ~isnan(d) && ~isnan(a)        % Checking that the angle and distance is acceptable
                Mn=mean(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Intn(d)=Intn(d)+Mn;
                Intd(d)=Intd(d)+Md;
                V_normal(d)=V_normal(d)+Vn;
                V_dark(d)=V_dark(d)+Vd;
                V_shot(d)=V_shot(d)+Md;
                V_corr(d)=V_corr(d)+(Vd-Mdv-Md);

                Compteur(d)=Compteur(d)+1;
            end
        end
    end
end

for d=1:Nmb_of_points
    Intn(d)=Intn(d)/Compteur(d);
    Intd(d)=Intd(d)/Compteur(d);
    V_normal(d)=V_normal(d)/Compteur(d);
    V_dark(d)=V_dark(d)/Compteur(d);
    V_shot(d)=V_shot(d)/Compteur(d);
    V_corr(d)=V_corr(d)/Compteur(d);

    Speckle(d)=V_normal(d)/Intn(d)^2;
    Speckle_corr(d)=V_corr(d)/Intd(d)^2;
end

%% Plots

t=tiledlayout(2,2);
title(t,'Computations for speckle contrast')

nexttile
plot(Distances_plot,Intd,"+")
xlabel('S-D separation (mm)')
ylabel('Intensity (e-)')
nexttile
plot(Distances_plot,log(Mean_dark_V),'+',Distances_plot,log(V_dark),"*",Distances_plot,log(V_shot),"x",Distances_plot,log(V_corr),"o")
xlabel('S-D separation (mm)')
ylabel('log(Variance) (e-^2)')
legend("Dark Variance","Variance of dark substracted intensity","Shot noise variance","Corrected variance")
nexttile
plot(Intd,log(Mean_dark_V),'+',Intd,log(V_dark),"*",Intd,log(V_shot),"x",Intd,log(V_corr),"o")
xlabel('Intensity (e-)')
ylabel('log(Variance) (e-^2)')
legend("Dark Variance","Variance of dark substracted intensity","Shot noise variance","Corrected variance")
nexttile
plot(Distances_plot,Speckle,"+",Distances_plot,Speckle_corr,"*")
xlabel('S-D separation (mm)')
ylabel('K^2')
legend("Uncorrected speckle contrast","Corrected speckle contrast")