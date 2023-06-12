%%
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec02.mp4');
Source1=mean(frames,4);
Source1=imresize(Source1,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec03.mp4');
Source2=mean(frames,4);
Source2=imresize(Source2,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec04.mp4');
Source3=mean(frames,4);
Source3=imresize(Source3,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec05.mp4');
Source4=mean(frames,4);
Source4=imresize(Source4,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec06.mp4');
Source5=mean(frames,4);
Source5=imresize(Source5,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec07.mp4');
Source6=mean(frames,4);
Source6=imresize(Source6,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec08.mp4');
Source7=mean(frames,4);
Source7=imresize(Source7,[350 700]);
frames=acquisition('C:\Users\qviot\Documents\MATLAB\Real_Data\230525\rec09.mp4');
Source8=mean(frames,4);
Source8=imresize(Source8,[350 700]);


%%
Source=Source1+Source2+Source3+Source4+Source5+Source6+Source7+Source8;
Source=mean(Source,4);
Source=Source/max(Source,[],'all');
% Source=Source*4096;
image(Source);

%%
x=linspace(0,300*(3.45*10^(-3)),301);
X=improfile(Source,[192,492],[171,171]);
X=squeeze(X);
y=linspace(0,300*(3.45*10^(-3)),301);
Y=improfile(Source,[342,342],[20,320]);
Y=squeeze(Y);
plot(x,X,y,Y)
xlabel('Distance in mm')
ylabel('Intensity')
legend("Horizontal profile","Vertical profile")