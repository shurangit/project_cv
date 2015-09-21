close all;
clc;
%read images and make them grayscale;
n=2;
sigma=1;
w=ceil(2.5*sigma)*2+1;
for i=1:n
str=strcat('DSC_028 (',int2str(i),').jpg');
str=strcat('DanaHallWay1\',str);
Io(:,:,:,i)=imread(str);
I(:,:,i)=rgb2gray(Io(:,:,:,i));
%compute Sx, Sy%
filter=fspecial('gaussian',[w,w], sigma);
filter_y=[-1 0 1];
filter_x=[-1 0 1]';
I_g=imfilter(double(I(:,:,i)),filter,'replicate');
[M,N]=size(I_g);
S_x=imfilter(I_g,filter_x,'replicate');
S_y=imfilter(I_g,filter_y,'replicate');
%define M matrix%
S2x=S_x.^2;
S2y=S_y.^2;
Sxy=S_x.*S_y;
w_s1=7;
filter_s=1/(w_s1*w_s1)*ones(w_s1,w_s1);
m=zeros(2,2,M,N);
m(1,1,:,:)=imfilter(S2x,filter_s,'replicate');
m(1,2,:,:)=imfilter(Sxy,filter_s,'replicate');
m(2,1,:,:)=m(1,2,:,:);
m(2,2,:,:)=imfilter(S2y,filter_s,'replicate');
%Compute R of each pixel%
R=zeros(M,N);
k=0.05;
for j=1:M
for jj=1:N
R(j,jj)=det(m(:,:,j,jj))-trace(m(:,:,j,jj))^2*k;
end
end
T_R=1000;
% Non-max suppression%
R1=R;
R1(R1<T_R)=0;
f(:,:,1)=[-1 0 0;0 1 0;0 0 0];
f(:,:,2)=[0 -1 0;0 1 0;0 0 0];
f(:,:,3)=[0 0 -1;0 1 0;0 0 0];
f(:,:,4)=[0 0 0;-1 1 0;0 0 0];
f(:,:,5)=[0 0 0;0 1 -1;0 0 0];
f(:,:,6)=[0 0 0;0 1 0;-1 0 0];
f(:,:,7)=[0 0 0;0 1 0;0 -1 0];
f(:,:,8)=[0 0 0;0 1 0;0 0 -1];
for iii=1:8
R2(:,:,iii)=imfilter(R1,f(:,:,iii));
end
for iii=1:8
R1(R2(:,:,iii)<0)=0;
end
I_R1=zeros(M,N);
%I_R1(R>T_R)=1;
I_R1(R1>0)=1;
%imtool images and corner is 1%
I_R(:,:,i)=I_R1;
II=I(:,:,i);
II(I_R1==1)=0;
imtool(II)
end
x=cell(2,1);
y=cell(2,1);
for i=1:n
[xx, yy]=find(I_R(:,:,i)==1);
x{i}=xx;
y{i}=yy;
end
w_s2=25;
w_r=(w_s2-1)/2;
p1=zeros(w_s2,w_s2);
p2=zeros(w_s2,w_s2);
[L1,L11]=size(x{1});
[L2,L22]=size(x{2});
sigmaF=1;
w_s3=ceil(2.5*sigmaF)*2+1;
filter_F=fspecial('gaussian',[w_s3,w_s3], sigmaF);
I_F(:,:,1)=imfilter(double(I(:,:,1)),filter_F,'replicate');
I_F(:,:,2)=imfilter(double(I(:,:,2)),filter_F,'replicate');
x1=x{1};
x2=x{2};
y1=y{1};
y2=y{2};
TR_2=0.995; %NCC threshold value
s=zeros(L1,L2);
a1=zeros(L1,6);
%Compute NCC%
for i=1:L1
p1=zeros(w_s2,w_s2);
if (x1(i)>w_r && x1(i)+w_r<M && y1(i)>w_r && y1(i)+w_r<N)
p1= I_F(x1(i)-w_r:x1(i)+w_r,y1(i)-w_r:y1(i)+w_r,1);
end
for j=1:L2
p2=zeros(w_s2,w_s2);
if (x2(j)>w_r && x2(j)+w_r<M && y2(j)>w_r && y2(j)+w_r<N)
p2= I_F(x2(j)-w_r:x2(j)+w_r,y2(j)-w_r:y2(j)+w_r,2);
end
p1v=p1(:);
p2v=p2(:);
if(p1v==0)
p1n=p1v;
else
p1n=p1v/norm(p1v);
end
if(p2v==0)
p2n=p2v;
else
p2n=p2v/norm(p2v);
end
p=p1n.*p2n;
s(i,j)=sum(p);
end
a=find(s(i,:)==max(s(i,:)));
a1(i,1)=a(1);
a1(i,2)=max(s(i,:));
a1(i,3)=x1(i);
a1(i,4)=y1(i);
a1(i,5)=x2(a(1));
a1(i,6)=y2(a(1));
if a1(i,2)<TR_2
a1(i,:)=0;
end
end
%combine 2 original images%
I1=[Io(:,:,:,1) Io(:,:,:,2)];
s1=size(a1);
imshow(I1);
for i=1:s1(1)
if a1(i,1)~=0
line([a1(i,4),a1(i,6)+N],[a1(i,3),a1(i,5)]);
end
end
c=a1((a1(:,1)~=0),3:6);
[sizea, sizea1]=size(c);
% Homography %
o_r=0.1;
Th_d=1.5;
i=5000;
while(i==5000)
for i=1:min(nchoosek(sizea,4),5000)
cp=zeros(4,4);
rd=[0 0 0 0];
while (rd(1)==rd(2) || rd(1)==rd(3) || rd(1)==rd(4) || rd(2)==rd(3) || rd(2)==rd(4) || rd(3)==rd(4))
rd=round((sizea-1)*rand(1,4)+1);
end
for j=1:4
cp(j,:)=c(rd(j),:);
end
cp_d=(cp(1,3:4)-cp(2,3:4)).*(cp(2,3:4)-cp(4,3:4)).*(cp(1,3:4)-cp(3,3:4)).*(cp(1,3:4)-cp(4,3:4)).*(cp(2,3:4)-cp(3,3:4)).*(cp(3,3:4)-cp(4,3:4));
if (cp_d(1)==0 && cp_d(2)==0)
continue
end
A=[cp(1,1) cp(1,2) 1 0 0 0 -cp(1,1)*cp(1,3) -cp(1,2)*cp(1,3);
0 0 0 cp(1,1) cp(1,2) 1 -cp(1,1)*cp(1,4) -cp(1,2)*cp(1,4);
cp(2,1) cp(2,2) 1 0 0 0 -cp(2,1)*cp(2,3) -cp(2,2)*cp(2,3);
0 0 0 cp(2,1) cp(2,2) 1 -cp(2,1)*cp(2,4) -cp(2,2)*cp(2,4);
cp(3,1) cp(3,2) 1 0 0 0 -cp(3,1)*cp(3,3) -cp(3,2)*cp(3,3);
0 0 0 cp(3,1) cp(3,2) 1 -cp(3,1)*cp(3,4) -cp(3,2)*cp(3,4);
cp(4,1) cp(4,2) 1 0 0 0 -cp(4,1)*cp(4,3) -cp(4,2)*cp(4,3);
0 0 0 cp(4,1) cp(4,2) 1 -cp(4,1)*cp(4,4) -cp(4,2)*cp(4,4)];
b=[cp(1,3) cp(1,4) cp(2,3) cp(2,4) cp(3,3) cp(3,4) cp(4,3) cp(4,4)]';
if (det(A)<10)
continue
end
h=A\b;
H=[h(1:3)';h(4:6)';h(7:8)',1];
num_out=0;
for j=1:sizea
co1=[c(j,1),c(j,2),1]';
co2(:,j)=H*co1;
co2(1,:)=co2(1,:)./(co2(3,:));
co2(2,:)=co2(2,:)./(co2(3,:));
co2(3,:)=co2(3,:)./(co2(3,:));
d(j)=((co2(1,j)-c(j,3))^2+(co2(2,j)-c(j,4))^2)^0.5;
if d(j)>Th_d
num_out=num_out+1;
end
end
out_r=num_out/sizea;
if out_r<o_r
i
out_r
break
end
i
end
o_r=o_r+0.05;
Th_d=1.5*Th_d;
o_r
Th_d
end
%Wrap 2 images%
I_out=zeros(2*M,2*N,3);
I_out(0.5*M+1:1.5*M,0.5*N+1:1.5*N,:)=double(Io(:,:,:,1));
for i=1:2*M
for j=1:2*N
Pc1=[i-0.5*M j-0.5*N 1]';
Pc2=H*Pc1;
Pc2=Pc2./Pc2(3);
cx=Pc2(1);
cy=Pc2(2);
if (cx>0.5 && cx<M+0.5 && cy>0.5 && cy<N+0.5)
cx2=ceil(cx);
cx1=cx2-1;
if cx<1
cx2=1;
cx1=1;
end
if cx>M
cx2=M;
cx1=M;
end
cy2=ceil(cy);
cy1=cy2-1;
if cy<1
cy2=1;
cy1=1;
end
if cy>N
cy2=N;
cy1=N;
end
dx2=abs(cx2-cx);
dx1=1-dx2;
dy2=abs(cy2-cy);
dy1=1-dy2;
PI2_11=double(Io(cx1,cy1,:,2));
PI2_12=double(Io(cx1,cy2,:,2));
PI2_21=double(Io(cx2,cy1,:,2));
PI2_22=double(Io(cx2,cy2,:,2));
PI2=dy2*(dx2*PI2_11+dx1*PI2_21)+dy1*(dx2*PI2_12+dx1*PI2_22);
w2=min([cx-1,M-cx,cy-1,N-cy]);
PI1=zeros(1,1,3);
w1=0;
if (i>0.5*M && i<1.5*M+1 && j>0.5*N && j<1.5*N+1)
PI1=I_out(i,j,:);
w1=min([i-0.5*M-1,1.5*M-i,j-0.5*N-1,1.5*N-j]);
end
I_out(i,j,:)=(w1^2.*PI1+w2^2.*PI2)/(w1^2+w2^2);
end
end
end
row1=0.5*M;
row2=0.5*M;
col1=0.5*N;
col2=0.5*N;
for i=1:0.5*M
if I_out(i,:,:)==0
row1=row1-1;
end
end
for i=1.5*M+1:2*M
if I_out(i,:,:)==0
row2=row2-1;
end
end
for i=1:0.5*N
if I_out(:,i,:)==0
col1=col1-1;
end
end
for i=1.5*N+1:2*N
if I_out(:,i,:)==0
col2=col2-1;
end
end
I_out=uint8(I_out(0.5*M-row1+1:1.5*M+row2,0.5*N-col1+1:1.5*N+col2,:));
imtool(I_out);