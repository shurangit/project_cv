%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diparity map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
imtool close all
close all

sigma=1;
w=ceil(2.5*sigma)*2+1;
n=2;
fn=cell(1,n);
Io=cell(1,n);
I=cell(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read images
uiwait(msgbox('Please select two images for Mosaicing','Image Mosaicing'));
for ii=1:n;
[fname,pname] = uigetfile('*.jpg','Select an image');

if fname==0
   IM=0;
   return 
end

fn{ii}=strcat(pname,fname);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Harris Corner 
for i=1:n
     Io{i}=imread(fn{i}); 
      %%%%%%%%%%%%%%%
         %if i==1
      %    Io{1}=Io{1}(50:250,50:350,:);
         %end
      %%%%%%%%%%%%%%
     I{i}=rgb2gray(Io{i}); 
  
   filter=fspecial('gaussian',[w,w], sigma);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gaussian filter

   filter_Gy=[-1 0 1];
   filter_Gx=[-1 0 1]';
   I_f=imfilter(double(I{i}),filter,'replicate');
   [M(i),N(i)]=size(I_f);
   I_x=imfilter(I_f,filter_Gx,'replicate');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient
   I_y=imfilter(I_f,filter_Gy,'replicate');
   
   I2x=I_x.^2;
   I2y=I_y.^2;
   Ixy=I_x.*I_y;
   
   w_s=7;
   filter_s=1/(w_s*w_s)*ones(w_s,w_s);
   m=zeros(2,2,M(i),N(i));
   m(1,1,:,:)=imfilter(I2x,filter_s,'replicate');
   m(1,2,:,:)=imfilter(Ixy,filter_s,'replicate');
   m(2,1,:,:)=m(1,2,:,:);
   m(2,2,:,:)=imfilter(I2y,filter_s,'replicate');
     
  R=zeros(M(i),N(i));
  k=0.05;
   for j=1:M(i)
       for jj=1:N(i)
   R(j,jj)=det(m(:,:,j,jj))-trace(m(:,:,j,jj))^2*k;
   
       end
   end
   
   T_R=10000;
   %%%%%%%%%%%%%%%%%%%%%%%%%%5
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
   
   R2=zeros(M(i),N(i),8);
   for iii=1:8
      R2(:,:,iii)=imfilter(R1,f(:,:,iii)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonmax suppression
   end
   for iii=1:8
      R1(R2(:,:,iii)<0)=0;
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%  
   I_R1=zeros(M(i),N(i));
   I_R1(R1>0)=1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   I_R{i}=I_R1;
  II=I{i};
  II(I_R1==1)=0;
   imshow(Io{i})%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot corners
   [X,Y]=find(I_R1==1);
   hold on
   plot(Y,X,'.')
   figure
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matching 
x=cell(2,1);
y=cell(2,1);
for i=1:n
[xx, yy]=find(I_R{i}==1);
x{i}=xx;
y{i}=yy;
end

w_m=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_r=(w_m-1)/2;
p1=zeros(w_m,w_m);
p2=zeros(w_m,w_m);
[L1,L11]=size(x{1});
[L2,L22]=size(x{2});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaF=1;
w_F=ceil(2.5*sigmaF)*2+1;
 filter_F=fspecial('gaussian',[w_F,w_F], sigmaF);
I_F{1}=imfilter(double(I{1}),filter_F,'replicate');
I_F{2}=imfilter(double(I{2}),filter_F,'replicate');

x1=x{1};
x2=x{2};
y1=y{1};
y2=y{2};
%%%%%%%%%%%%%%%%%%%%%%
T2=0.995;
%%%%%%%%%%%%%%%%%%%%%%
s=zeros(L1,L2);
a1=zeros(L1,6);
for i=1:L1
    i
    L1
    p1=zeros(w_m,w_m);
      if (x1(i)>w_r &&  x1(i)+w_r<M(1) && y1(i)>w_r &&  y1(i)+w_r<N(1))
        p1= I_F{1}(x1(i)-w_r:x1(i)+w_r,y1(i)-w_r:y1(i)+w_r);
      
      end
    
      for j=1:L2
         p2=zeros(w_m,w_m);
         if (x2(j)>w_r &&  x2(j)+w_r<M(2) && y2(j)>w_r &&  y2(j)+w_r<N(2))
              p2= I_F{2}(x2(j)-w_r:x2(j)+w_r,y2(j)-w_r:y2(j)+w_r);
         end
         
         p1v=p1(:);
         p2v=p2(:);
         if(p1v==0)
            p1n=p1v; 
         else
         p1n=p1v/norm(p1v);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization
         end
        if(p2v==0)
            p2n=p2v; 
         else
         p2n=p2v/norm(p2v);
         end
         p=p1n.*p2n;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NCC
         s(i,j)=sum(p);
      end
      
      
     a=find(s(i,:)==max(s(i,:)));
     a1(i,1)=a(1);
     a1(i,2)=max(s(i,:));
     a1(i,3)=x1(i);
     a1(i,4)=y1(i);
     a1(i,5)=x2(a(1));
     a1(i,6)=y2(a(1));
     
     if a1(i,2)<T2
      a1(i,:)=0;
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I1=zeros(max(M),N(1)+N(2),3,'uint8');
for i=1:max(M)
  if i<=min(M)
    I1(i,:,:)=[Io{1}(i,:,:),Io{2}(i,:,:)];
  else
      if M(1)>M(2)
        I1(i,:,:)=[Io{1}(i,:,:),255*ones(1,N(2),3,'uint8')];
      else
        I1(i,:,:)=[255*ones(1,N(1),3,'uint8'),Io{2}(i,:,:)];   
      end
  end
end

s1=size(a1);

imshow(I1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show matching
for i=1:s1(1)
    if a1(i,1)~=0
        hold on
      plot([a1(i,4),a1(i,6)+N(1)],[a1(i,3),a1(i,5)],'r');    
    end
end

c=a1((a1(:,1)~=0),3:6);

[sizea, sizea1]=size(c);
if sizea<8
  uiwait(msgbox('error'));
  IM=0;
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_r=0.1;
Th_d=0.5;
i=min(nchoosek(sizea,8),5000);
while(i==min(nchoosek(sizea,8),5000))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANSAC
for i=1:min(nchoosek(sizea,8),5000) 
cp=zeros(8,4);

%rd=[0 0 0 0 0 0 0 0];
A=[];
c_in=zeros(sizea,sizea1);
%while (rd(1)==rd(2) || rd(1)==rd(3) || rd(1)==rd(4) || rd(2)==rd(3) || rd(2)==rd(4) || rd(3)==rd(4))
rd=round((sizea-1)*rand(1,8)+1);
%end

for j=1:8
   cp(j,:)=c(rd(j),:); 
end

for kk=1:8
A=[A;
    cp(kk,1)*cp(kk,3) cp(kk,1)*cp(kk,4) cp(kk,1) cp(kk,2)*cp(kk,3) cp(kk,2)*cp(kk,4) cp(kk,2) cp(kk,3) cp(kk,4) 1];
end
[U,S,V] = svd(A);
 F=V(:,end);
 f=[F(1:3),F(4:6),F(7:9)];
 [Uf,Sf,Vf] = svd(f);

Sf(3,3)=0;
f=Uf*Sf*Vf';

f1=f(1:2,:)';
f2=-f(3,:)';
er=f1\f2;
er=[er;1];
f3=f(:,1:2);
f4=-f(:,3);
el=f3\f4;
el=[el;1];

num_out=0;
d=zeros(sizea,1);
for j=1:sizea
    co1=[c(j,1),c(j,2),1]';
   
    lr=f*co1;
    d(j)=abs(lr(1)*c(j,3)+lr(2)*c(j,4)+lr(3))/(lr(1)^2+lr(2)^2)^0.5;
    if d(j)>Th_d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       num_out=num_out+1;
    else
        c_in(j,:)=c(j,:);
    end
end

out_r=num_out/sizea;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if out_r<o_r
    i
    out_r
    break
end
i
end
o_r
Th_d
o_r=o_r+0.05;
Th_d=0+Th_d;

end

for i=1:(sizea-num_out)
        hold on
        plot([c_in(i,2),c_in(i,4)+N(1)],[c_in(i,1),c_in(i,3)]);    
end
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Best FMatrix
[M1,N1]=size(c_in);

A2=[];


for i=1:M1
A2=[A2;
    c_in(i,1)*c_in(i,3) c_in(i,1)*c_in(i,4) c_in(i,1) c_in(i,2)*c_in(i,3) c_in(i,2)*c_in(i,4) c_in(i,2) c_in(i,3) c_in(i,4) 1];
end
[U,S,V] = svd(A2);
 F=V(:,end);
 fff=[F(1:3),F(4:6),F(7:9)];
 [Uf,Sf,Vf] = svd(f);

Sf(3,3)=0;
fff=Uf*Sf*Vf';

f1=f(1:2,:)';
f2=-f(3,:)';
er=f1\f2;
er=[er;1];
f3=f(:,1:2);
f4=-f(:,3);
el=f3\f4;
el=[el;1];

in_s=1;
search_step=1;
Iddx=ones(ceil(M(1)/in_s)-1,ceil(N(1)/in_s)-1)*1000;
Iddy=ones(ceil(M(1)/in_s)-1,ceil(N(1)/in_s)-1)*1000;
w_r=(15-1)/2;
cos_window = hann(w_r*2+1) * hann(w_r*2+1)';
%cos_window =ones(w_r*2+1,w_r*2+1);
 
for i=2:ceil(M(1)/in_s)-1
   for j=2:ceil(N(1)/in_s)-1
    co1=[1+in_s*i,1+j*in_s,1]';
    if (co1(1)>w_r &&  co1(1)+w_r<M(1) && co1(2)>w_r &&  co1(2)+w_r<N(1))
        p1= I_F{1}(co1(1)-w_r:co1(1)+w_r,co1(2)-w_r:co1(2)+w_r);
        p1=p1.*cos_window;
    else
        continue
    end
    p1v=p1(:);
    if(p1v==0)
        p1n=p1v; 
    else
        p1n=p1v/norm(p1v);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization
    end
     lr=fff*co1;
     %[FF, epipolarInliers, status] = estimateFundamentalMatrix(...
  %c_in(:,1:2), c_in(:,3:4), 'Method', 'RANSAC', ...
  %'NumTrials', 10000, 'DistanceThreshold', 0.1, 'Confidence', 99.99);
    % lr=FF*co1;
     %f1=FF(1:2,:)';
     %f2=-FF(3,:)';
     %er=f1\f2;
    % er=[er;1];
    % f3=FF(:,1:2);
    % f4=-FF(:,3);
    % el=f3\f4;
     %el=[el;1];
   %  imshow(I1)
    % hold on
     if abs(er(3)/er(2))<0.001
    %   plot([N(1),N(2)+N(1)],[-lr(3)/lr(1),-lr(3)/lr(1)]);  
     else
     %plot([N(1),er(2)+N(1)],[-lr(3)/lr(2),er(1)]);
     end
   %  plot([co1(2),er(2)],[co1(1),er(1)]);
  
    co2=[0 0 1]';
    co22=co2;
    
   
    
   

          tty=co1(2)-80;
          ttty=co1(2)-10;
          ttx=round(-lr(2)/lr(1)*tty-lr(3)/lr(1));
          tttx=round(-lr(2)/lr(1)*ttty-lr(3)/lr(1));
          s0=0;
   for jj=max(1,tty):search_step:min(ttty,N(2))
           co2(2)=jj;
           co2(1)=round(-lr(2)/lr(1)*jj-lr(3)/lr(1));
         if (co2(1)>w_r &&  co2(1)+w_r<M(2) && co2(2)>w_r &&  co2(2)+w_r<N(2))
              p2= I_F{2}(co2(1)-w_r:co2(1)+w_r,co2(2)-w_r:co2(2)+w_r);
              p2=p2.*cos_window;
         else
             continue
         end
          p2v=p2(:);
         if(p2v==0)
            p2n=p2v; 
         else
         p2n=p2v/norm(p2v);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization
         end
     p=p1n.*p2n;    
     ss=sum(p);
    
     i
     j
     if ss>s0 && ss>0.95
          % plot(co1(2),co1(1),'+');
          % plot(co2(2)+N(1),co2(1),'+');
          
           if j>2
             if (co2(2)-co1(2))-(Iddy(i,j-1))>15 && abs(Iddy(i,j-1)-Iddy(i,j-2))<15 &&  Iddy(i,j-1)~=1000  &&  Iddy(i,j-2)~=1000
                break 
             end
           end
         co22=co2;
        Iddx(i,j)=co2(1)-co1(1);
        Iddy(i,j)=co2(2)-co1(2);
        s0=ss;
     end   
         
   end
  
   
   
  
 % plot(co22(2)+N(1),co22(1),'o');
 % plot(co1(2),co1(1),'o');
%  p3=[p1;p2];
%  imshow(p3/255)
  s0
  if i>34 && i<38
      
  end
  
   
  % drawnow
   end
end

[xi, yi] = meshgrid(1:N(1), 1:M(1));
Iddy2=interp2(Iddy, xi/in_s, yi/in_s);
Iddx2=interp2(Iddx, xi/in_s, yi/in_s);