
clear all
close all
base_path = './data/';

sigma=1;
sigma_c=2;
patch_wise=[21,21];
K=10;               %num of clusters
thrd_sz=3;
thrd=3;
target_sz=[40,100];

[word_dscrpt,word_dspl,word_sumd]=OCD_Train(base_path,patch_wise,K,sigma_c,sigma);
% saving training data
 save (strcat(base_path,'dscrpt.mat'),'word_dscrpt');
 save (strcat(base_path,'dspl.mat'),'word_dspl');
 save (strcat(base_path,'sumd.mat'),'word_sumd');

