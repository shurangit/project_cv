function [ c_pos] = OCD_test( im,word_dscrpt,word_dspl, word_sumd,patch_wise,sigma_c,sigma,thrd_sz,thrd,target_sz )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

w_s=7;
T_R=1000;

patch_half=floor(patch_wise/2);
im2=im;
if size(im,3) > 1,
		im = rgb2gray(im);
end
%BG=zeros(600,600);

im=double(im);
[M,N]=size(im);
pos = Harris_Corner(im,sigma_c,w_s,T_R);
   
    X=pos(:,1);
    Y=pos(:,2);
      %imshow(im/255)
    %hold on
 %  plot(Y,X,'+')
    %pnts{i}=zeros(numel(X),patch_wise^2);
   % pnts_info{i}=zeros(numel(X),4);
   
    w=ceil(2.5*sigma)*2+1;
   filter=fspecial('gaussian',[w,w], sigma);
   
   im=imfilter(double(im),filter,'replicate'); 
      %  close all
     %imshow(im2);
   pnts_info=[];
   dst=zeros(1,size(word_dspl,1));
   b=1:size(word_dspl,1);
    for j=1:numel(X)
       try 
           a= im(X(j)-patch_half(1):X(j)+patch_half(1),Y(j)-patch_half(2):Y(j)+patch_half(2));
           a=double(a(:)');
           if(sum(a)~=0)
            a=a/norm(a);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization
           else
               a
           end
       for i=1:size(word_dspl,1)
          dst(i)=dist(word_dscrpt(i,:),a');
       end
       min_dist=min(dst);
       dd=dst(dst>min_dist);
       secd_mindist=min(dd);
       if secd_mindist>1.1*min_dist || secd_mindist<0.3
           secd_mindist=0;
       end
       dspl=word_dspl(dst==min_dist,:);
       n=b(dst==min_dist);
       n=n(1);
       dspl=dspl(1,:);
       s_d= word_sumd(dst==min_dist);
       s_d=s_d(1);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if secd_mindist~=0
         dspl2=word_dspl(dst==secd_mindist,:);
     
       n2=b(dst==secd_mindist);
       n2=n2(1);
       dspl2=dspl2(1,:);
       s_d2= word_sumd(dst==secd_mindist);
       s_d2=s_d2(1);  
       end
       %s=[s;a];    
       %pnts{i}(j,:)=a;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if min_dist<0.5
         min_dist=1;
         pnts_info=[pnts_info;[X(j),Y(j)]+dspl([1,2]), n,1/(s_d*min_dist)];
         pnts_info=[pnts_info;[X(j),Y(j)]+dspl([3,4]), n,1/(s_d*min_dist)];
        %pnts_info=[pnts_info;[X(j),Y(j)]+dspl([1,2]), n,100-s_d];
        %pnts_info=[pnts_info;[X(j),Y(j)]+dspl([3,4]), n,100-s_d];
         rect_position=[[Y(j),X(j)]-patch_wise([2,1])/2,patch_wise([2,1])];
         if secd_mindist>0
         pnts_info=[pnts_info;[X(j),Y(j)]+dspl2([1,2]), n2,1/(s_d2*min_dist)];
         pnts_info=[pnts_info;[X(j),Y(j)]+dspl2([3,4]), n2,1/(s_d2*min_dist)];      
         end
         
         if  n2==1
         
       %  rectangle('Position',rect_position, 'EdgeColor','r'); 
      %  hold on
       %  plot([Y(j),Y(j)+dspl2(2)],[X(j),X(j)+dspl2(1)],'LineWidth',2);
       %  plot([Y(j),Y(j)+dspl2(4)],[X(j),X(j)+dspl2(3)],'lineWidth',2);
       %  hold on
       %  n
         
        %s_d
       % min_dist
        
        %1/(s_d*min_dist)
      
           
        end
        
        
       end
       catch
           
           continue
       end  
    end

b_img=zeros(M,N);
for i=1:size(pnts_info,1)
    x=round(pnts_info(i,1));
    y=round(pnts_info(i,2));
   if x>0 && x<M+1 && y>0 && y<N+1
    b_img(x,y)=b_img(x,y)+pnts_info(i,4);
      %b_img(x,y)=b_img(x,y)+1;
   end
end
filter=ones(thrd_sz);
responses=imfilter(b_img,filter);

 f(:,:,1)=[-1 0 0;0 1 0;0 0 0];
   f(:,:,2)=[0 -1 0;0 1 0;0 0 0];
   f(:,:,3)=[0 0 -1;0 1 0;0 0 0];
   f(:,:,4)=[0 0 0;-1 1 0;0 0 0];
   f(:,:,5)=[0 0 0;0 1 -1;0 0 0];
   f(:,:,6)=[0 0 0;0 1 0;-1 0 0];
   f(:,:,7)=[0 0 0;0 1 0;0 -1 0];
   f(:,:,8)=[0 0 0;0 1 0;0 0 -1];
   
   r2=zeros(M,N,8);
   for iii=1:8
      r2(:,:,iii)=imfilter(responses,f(:,:,iii)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonmax suppression
   end
   for iii=1:8
      responses(r2(:,:,iii)<0)=0;
   end
  % close
  % imshow(BG);
   imshow(im2)
   pos_target=[];
for i=1:M
    for j=1:N
     if responses(i,j)>thrd
     %rect_position=[[j,i]-target_sz([2,1])/2,target_sz([2,1])];
     %rectangle('Position',rect_position, 'EdgeColor','r'); 
     pos_target=[pos_target; i,j];
    
     end
    end
end  
      c_pos=[];
      
     if size(pos_target,1)>2
     [cidx, c_pos] = kmeans(pos_target, 3,'EmptyAction','drop');
     c_pos(isnan(c_pos)) = [];
     c_pos=[c_pos([1,3]);c_pos([2,4])];
     elseif size(pos_target,1)>1
     [cidx, c_pos] = kmeans(pos_target, 2,'EmptyAction','drop');
     c_pos(isnan(c_pos)) = [];
     elseif size(pos_target,1)==1
         c_pos=pos_target;
     end
           if size(c_pos,2)>2
         
           end
      if size(c_pos,1)==3
         pa=c_pos(1,:);
         pb=c_pos(2,:);
         pc=c_pos(3,:);
         
          if dist(pa,pb')<target_sz(1)/3
           c_pos=[(pa+pb)/2;pc];
          elseif dist(pa,pc')<target_sz(1)/3
           c_pos=[(pa+pc)/2;pb];
          elseif dist(pb,pc')<target_sz(1)/3
           c_pos=[(pb+pc)/2;pa];
          end     
           if size(c_pos,2)>2
         
           end
     end
     if size(c_pos,1)==2
         pa=c_pos(1,:);
         pb=c_pos(2,:);
         if dist(pa,pb')<target_sz(1)/3
          c_pos=(pa+pb)/2;
         end
     end
   
     if size(c_pos,2)>2
         
     end
     
     for i=1:size(c_pos,1)
       rect_position=[[c_pos(i,2),c_pos(i,1)]-target_sz([2,1])/2,target_sz([2,1])];
       rectangle('Position',rect_position, 'EdgeColor','r','LineWidth',2);
     end
end

