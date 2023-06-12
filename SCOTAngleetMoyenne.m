%% Getting matrices from .mj2 or mp4. file

clear variables
Normal_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230609\rec05.mp4');
Dark_frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230609\rec06.mp4');
Normal_frames=double(Normal_frames(:,:,1:500));
Dark_frames=double(Dark_frames(:,:,1:100));
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

%% Computation of distances and K

% Counts in electrons
CF=2.51;
Normal_frames_e=Normal_frames*CF;
Dark_frames_e=Dark_frames*CF;

% Substracting by mean of dark frames
Dark_mean=mean(Dark_frames_e,3);
Darksubstracted_frames=Normal_frames_e-Dark_mean;

% Dividing all the possible distances by intervals
Nmb_of_points_d=40;
Dist_min=0.5;
Dist_max=2;
Ecart_d=(Dist_max-Dist_min)/(Nmb_of_points_d);
Distances_plot=linspace(Dist_min+Ecart_d/2,Dist_max-Ecart_d/2,Nmb_of_points_d);

% Dividing all the possible angles by intervals
Nmb_of_points_a=8;
Angle_min=-pi;
Angle_max=pi;
Ecart_a=(Angle_max-Angle_min)/Nmb_of_points_a;
Angles_plot=linspace(Angle_min-Ecart_a/2,Angle_max-Ecart_a/2,Nmb_of_points_a);

% Storing the values for every acquisition
V_normal=zeros(Nmb_of_points_a,Nmb_of_points_d);
Intn=zeros(Nmb_of_points_a,Nmb_of_points_d);
Intd=zeros(Nmb_of_points_a,Nmb_of_points_d);
V_dark=zeros(Nmb_of_points_a,Nmb_of_points_d);
V_shot=zeros(Nmb_of_points_a,Nmb_of_points_d);
V_corr=zeros(Nmb_of_points_a,Nmb_of_points_d);
Speckle=zeros(Nmb_of_points_a,Nmb_of_points_d);
Speckle_corr=zeros(Nmb_of_points_a,Nmb_of_points_d);

% Size of the sliding pixel window
W_size=7;
RWindow=(W_size-1)/2; %Distance in pixels from center to the edge of the window

% Computing dark variance 
Mdv=0;
Compteur=0;
Mean_dark_V=zeros(Nmb_of_points_a,Nmb_of_points_d);
for l=1:L_dark
    Dark_frame=Dark_frames_e(:,:,l);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            Mdv=Mdv+var(Dark_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
            Compteur=Compteur+1;
        end
    end
end
Mdv=Mdv/Compteur;
Mean_dark_V(:)=Mdv;

% One computation of distance/angle for each source
% We attirbute each distance value to an indice which is stored in Distances_sources (same for angle)
Distances_sources=zeros(N_sources,Height-6,Width-6);
Angles_sources=zeros(N_sources,Height-6,Width-6);
for l=1:N_sources
    for i=1:Height-6
        for j=1:Width-6
            D=pdist2([Sources(l,1)*0.00943,Sources(l,2)*0.012],[(j+3)*0.00943,(i+3)*0.012]); % Distances in mm
            if D<Distances_plot(1)-(Ecart_d/2) || D>Distances_plot(end)+(Ecart_d/2)
                Distances_sources(l,i,j)=nan;
            else
                for c=1:Nmb_of_points_d 
                    if D>=Distances_plot(c)-(Ecart_d/2) && D<Distances_plot(c)+(Ecart_d/2)
                        Distances_sources(l,i,j)=c;
                    end
                end
            end
            A=angle(complex((j+3)-Sources(l,1),((i+3)-Sources(l,2))));
            for c=1:Nmb_of_points_a-1 
                if A>=Angles_plot(c) && A<Angles_plot(c+1)
                    Angles_sources(l,i,j)=c;
                end
            end
            if Angles_sources(l,i,j)==0
                Angles_sources(l,i,j)=Nmb_of_points_a;
            end
        end
    end
end

% Computation of variances and K
Compteur=zeros(Nmb_of_points_a,Nmb_of_points_d);
for l=1:L_normal
    Normal_frame=Normal_frames_e(:,:,l);
    Darksubstracted_frame=Darksubstracted_frames(:,:,l);
    n=mod(l-1,N_sources)+1;    % Gives the number of the source n for the frame l
    Angle_source=Angles_sources(n,:,:);
    Distance_source=Distances_sources(n,:,:);
    for i=1:Height-(W_size-1)
        for j=1:Width-(W_size-1)
            a=Angle_source(1,i,j);
            d=Distance_source(1,i,j);
            if ~isnan(d)
                Mn=mean(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Md=mean(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),"all");
                Vn=var(Normal_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Vd=var(Darksubstracted_frame(i:i+(W_size-1),j:j+(W_size-1)),0,"all");
                Intn(a,d)=Intn(a,d)+Mn;
                Intd(a,d)=Intd(a,d)+Md;
                V_normal(a,d)=V_normal(a,d)+Vn;
                V_dark(a,d)=V_dark(a,d)+Vd;
                V_shot(a,d)=V_shot(a,d)+Md;
                V_corr(a,d)=V_corr(a,d)+(Vd-Mdv-Md);

                Compteur(a,d)=Compteur(a,d)+1;
            end
        end
    end
end

for a=1:Nmb_of_points_a
    for d=1:Nmb_of_points_d
    Intn(a,d)=Intn(a,d)/Compteur(a,d);
    Intd(a,d)=Intd(a,d)/Compteur(a,d);
    V_normal(a,d)=V_normal(a,d)/Compteur(a,d);
    V_dark(a,d)=V_dark(a,d)/Compteur(a,d);
    V_shot(a,d)=V_shot(a,d)/Compteur(a,d);
    V_corr(a,d)=V_corr(a,d)/Compteur(a,d);

    Speckle(a,d)=V_normal(a,d)/Intn(a,d)^2;
    Speckle_corr(a,d)=V_corr(a,d)/Intd(a,d)^2;
    end
end

%% Plots

title('Intensity in fucntion of S-D separation for different angles')
hold on

for a=1:Nmb_of_points_a
    plot(Distances_plot,log(Intd(a,:)))
end 
xlabel('S-D separation (mm)')
ylabel('Log(intensity) (e-)')
ylim([2 5])
legend('180º','225º','270º','315º','0º','45º','90º','135º')
