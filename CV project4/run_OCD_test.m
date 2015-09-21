close all
base_path = './data/';
sigma_c=2;
sigma=1;
patch_wise=[21,21];
%K=7;
thrd_sz=5;
%thrd=0.32;

%thrd=0.15;
target_sz=[40,100];

try 
a1=word_dscrpt;
a1=word_dspl;
a1=word_sumd;
catch
   try 
   f= load([base_path,'sumd.mat']);
   word_sumd=f.word_sumd;
   f= load([base_path,'dscrpt.mat']);
   word_dscrpt=f.word_dscrpt;
   f= load([base_path,'dspl.mat']);
   word_dspl=f.word_dspl;
   catch
       disp('No training data');
       return
   end
end

thrd=(sum(1./(word_sumd))/numel(word_sumd)+1/min(word_sumd))*1.45;
[video_path,img_files] = load_video(base_path);
n=numel(img_files);
positions = cell(1,n);
AA=zeros(200,300)+0.5;
imshow(AA)
hold
for i=1:n
  im=imread([video_path img_files{i}]);   
  imshow(AA)
  %imshow(im)
[ pos_target] = OCD_test( im,word_dscrpt,word_dspl, word_sumd,patch_wise,sigma_c,sigma,thrd_sz,thrd,target_sz );
 positions{i}=pos_target;
pause(0.1)
end

gth=load([video_path,'CarsGroundTruthBoundingBoxes.mat']);
ground_truth=gth.groundtruth;
[rate_p,rate_ms]=show_precision_OCD(positions,target_sz, ground_truth, video_path,video_path,img_files)