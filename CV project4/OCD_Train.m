function [word_dscrpt,word_dspl, word_sumd]=OCD_Train(base_path,patch_wise,K,sigma_c,sigma)
%clear all
%base_path = './data/';

%sigma=1;
w_s=7;
T_R=10000;
%patch_wise=[21,21];
%K=10;
[video_path,img_files] = load_video(base_path);
 patch_half=floor(patch_wise/2);
num_img_files=numel(img_files);
%num_img_files=50;

s=[];
%pnts=cell(num_img_files);
pnts_info=[];
center=zeros(num_img_files,2);


for i=1:num_img_files
   im=imread([video_path img_files{i}]);
    if size(im,3) > 1,
		im = rgb2gray(im);
    end
    im=double(im);
   % im2=im;
    [M,N]=size(im);
    center(i,:)=floor([M/2,N/2]);
    pos = Harris_Corner(im,sigma_c,w_s,T_R);
   
     w=ceil(2.5*sigma)*2+1;
   filter=fspecial('gaussian',[w,w], sigma);
   
   im=imfilter(im,filter,'replicate'); 
    % imshow(im/255)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot corners
    X=pos(:,1);
    Y=pos(:,2);
   % hold on
  % plot(Y,X,'+')
    
    
  
    %pnts{i}=zeros(numel(X),patch_wise^2);
   % pnts_info{i}=zeros(numel(X),4);
    for j=1:numel(X)
       try 
           a= im(X(j)-patch_half(1):X(j)+patch_half(1),Y(j)-patch_half(2):Y(j)+patch_half(2));
           a=double(a(:)');
           if(sum(a)~=0)
            a=a/norm(a);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization
           else
               a
           end
       s=[s;a];    
       %pnts{i}(j,:)=a;
       pnts_info=[pnts_info;X(j),Y(j),i,j,0,0 0];
           
       catch
           
           continue
       end  
    end
    % pnts{i}=pnts{i}(pmt_info{i}~=0);
     %pnts{i}=pnts{i}(pmt_info{i}~=0);
end

[cidx, ctrs,SUMD,D] = kmeans(s, K);
for j=1:15
%if ~isnan(ctrs(:,1))
for i=1:size(cidx,1)
   c=cidx(i);
   d=D(i,c);
   if d>0.15
      s(i,:)=0; 
   end 
end
try
[cidx, ctrs,SUMD,D] = kmeans(s, K);
catch
    break
end
%[cidx, ctrs,SUMD,D] = kmeans(s, K,'EmptyAction','drop');
%else
 %   break
%end
end
pnts_info(:,5)=cidx;
n_o=numel(SUMD(SUMD==0));
TR=sum(SUMD)/(K-n_o)*1.2;
cidx1=1:K;
cidx1(SUMD>TR)=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cidx1(SUMD==0)=0;
pnts_info0=cell(1,num_img_files);
j=1;
k=1;
ff=size(pnts_info,1);
s_pos=cell(K,1);
num_pos=zeros(K,1)+eps;
for i=1:ff
   if   i==ff || pnts_info(min(i+1,ff),3)>k 
      p= pnts_info(j:i,:);
      d=D(j:i,:);
      pnts_info0{k}=[];
      pos0=center(k,:);
      for ii=1:K
      
       if cidx1(ii)>0
        try
        m= d(p(:,5)==ii,:);
          mm_min=min(m(:,ii));
          mm=m(m(:,ii)<min(1.1*mm_min,1),:);
          for iii=1:size(mm,1)
             num=p(d(:,ii)==mm(iii,ii),4);
             p1=p(p(:,4)==num,:);
             pos1=p1(:,1:2);
             p1(:,6:7)=pos0-pos1; 
            % s_pos(ii,:)=[s_pos(ii,:);pos0-pos1];
             s_pos{ii}=[s_pos{ii};pos0-pos1];
             num_pos(ii)=num_pos(ii)+1;
             pnts_info0{k}=[pnts_info0{k};p1(:,[1,2,5,6,7])];
          
          end
        catch
           continue 
        end
       end
      end
      j=i;
      k=k+1;
     % if k==num_img_files
        %  k=k-1;
     % end
          
   end
   
   
end
%avrg_pos=[s_pos(:,1)./(num_pos),s_pos(:,2)./(num_pos)];
avrg_pos=zeros(K,4);
for i=1:K
    if isempty(s_pos{i})
         avrg_pos(i,:)=0;
    else   
        [~, ct] = kmeans( s_pos{i}, 2);
        avrg_pos(i,1:2)=ct(1,:);
        avrg_pos(i,3:4)=ct(2,:);
    end
end

n=numel(avrg_pos(avrg_pos(:,1)~=0));
word_dscrpt=zeros(n,patch_wise(1)*patch_wise(2));
word_dspl=zeros(n,4);
word_sumd=zeros(n,1);
j=1;
for i=1:K
    if avrg_pos(i,1)~=0
        word_dscrpt(j,:)=ctrs(i,:);
        word_dspl(j,:)=avrg_pos(i,:);
        word_sumd(j)=SUMD(i);
        j=j+1;
    end
    
end
%end
%s=[];
%ctrs=[];
%pnts_info=[];
%pnts_info0=[];
%D=[];
close all
%for i=1:num_img_files
 %   n=size(pnts_info0{i},1);
  %  p=pnts_info0{i};
   %  im=imread([video_path img_files{i}]);
    % imshow(im);
    %for j=1:n
    %    pos=p(j,[1,2]);
    %    rect_position=[[pos(2),pos(1)]-patch_wise([2,1])/2,patch_wise([2,1])];
    %   
    %    rectangle('Position',rect_position, 'EdgeColor','g');
    %    p(j,3)
    %end
    
end