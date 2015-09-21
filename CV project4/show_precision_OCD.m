function [rp,rms]=show_precision_OCD(positions,target_sz, ground_truth, title,video_path,img_files)


	
	max_threshold = 1;  
	
	R=zeros(1000,1000);
    for i=1:1000
         for j=1:1000
          R(i,j)=(i-1)*1000+j;      
         end
    end
	if size(positions,2) ~= size(ground_truth,2),
		disp('Could not plot precisions, because the number of ground')
		disp('truth frames does not match the number of tracked frames.')
		return
    end
    n=0;
    nn=0;
    l=size(positions,2);
    for i=1:l
       n1=size(positions{i},1);
       p2=ground_truth(i).topLeftLocs;
       n2=size(p2,1);
       nn=nn+n1;
       im=imread([video_path img_files{i}]);   
       imshow(im);
       if i==99
           
       end
      for j=1:n1
           pos1=positions{i}(j,:);
           pos1=floor(pos1-target_sz/2);
            rect_position=[[pos1(2),pos1(1)],target_sz([2,1])];
       rectangle('Position',rect_position, 'EdgeColor','g','LineWidth',2);
       for jj=1:n2
           pos2=p2(jj,[2,1]);
            rect_position=[[pos2(2),pos2(1)],target_sz([2,1])];
       rectangle('Position',rect_position, 'EdgeColor','r','LineWidth',2);
           pos1
           pos2
           i
           if abs(pos2(1)-pos1(1))>target_sz(1)  ||   abs(pos2(2)-pos1(2))>target_sz(2)
               r=0;
               r
           else
               sz=target_sz-[abs(pos2(1)-pos1(1)) abs(pos2(2)-pos1(2))];
               a1=sz(1)*sz(2);
               a2=target_sz(1)*target_sz(2)*2-a1;
               r=a1/a2;
           end
            if r>0.5;
                 n=n+1; 
                 break
            end 
       end
      end
    end   
    rp=n/nn;
    
    
     n=0;
    nn=0;
    l=size(positions,2);
    for i=1:l
       n1=size(positions{i},1);
       p2=ground_truth(i).topLeftLocs;
       n2=size(p2,1);
       nn=nn+n2;
      for j=1:n2
           pos2=p2(j,[2,1]);
       for jj=1:n1
           pos1=positions{i}(jj,:);
           pos1=floor(pos1-target_sz/2);
           if abs(pos2(1)-pos1(1))>target_sz(1)  ||   abs(pos2(2)-pos1(2))>target_sz(2)
               r=0;
           else
               sz=target_sz-[abs(pos2(1)-pos1(1)) abs(pos2(2)-pos1(2))];
               a1=sz(1)*sz(2);
               a2=target_sz(1)*target_sz(2)*2-a1;
               r=a1/a2;
           end
            if r>0.5;
                 n=n+1; 
                 break
            end 
       end
      end
    end   
    rms=1-n/nn;
end