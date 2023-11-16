clear
close all
%path data
% path='/Users/platel/Desktop/exp/aurelie/to analyse/444118/clean_plane2/';
% load('/Users/platel/Desktop/exp/allnwb.mat')
% filename=cell2mat(nwbFiles (24));
% [path,name,ext] = fileparts(filename);
% path =[path '\'];
% path='/Users/platel/Desktop/exp/aurelie/to analyse/444113/444113_221010_plane0/';
% path='/Volumes/CossartNAS/Aurelie/brainbow analysis/444118/220912_plane0_000/';
% path='/Volumes/10.51.106.5/Data/Aurelie/brainbow analysis/444118/220912_plane1_000/';
path='/Volumes/10.51.106.5/Data/Aurelie/brainbow analysis/444175/221122_plane1/';

%path='/Users/platel/Desktop/dossier sans titre/';
% name='P2M_444113_221010_plane0_2023_03_06.12-20-10.nwb';
% name='P2M_444118_220912_000_plane0_2023_05_02.15-59-39.nwb';
% name='test';

%smooth image before ???
blue=imread([path,'blue.tif']);
red1040=imread([path,'red.tif']);
red890=imread([path,'red890.tif']);
green=imread([path,'green.tif']);
% calcium=imread([path,'calcium.tif']);
%composite=imread([path,'composite.tif']);

% PathSave='/Users/platel/Desktop/exp/analysis/';
% daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
% namefull=[PathSave daytime  '/'];
% mkdir (namefull)    % make folder for saving analysis
% path=namefull;

%import calcium image from Fall.mat suite2p
load([path 'Fall.mat'])
calcium=uint16(ops.refImg);
imwrite(calcium, ([path,'calcium.tif']), 'tif');
% calcium=ops.meanImgE;


%filter images (background substraction rolling ball): 
radius = 50;
se = strel('disk', radius);
%% 

% % Perform the rolling ball background subtraction
green = imtophat(green, se);
% green = imadjust(green)/4;
red1040 = imtophat(red1040, se);
% red1040 = imadjust(red1040)/2;
blue = imtophat(blue, se);
% blue =  imadjust(blue);
% calcium=imadjust(calcium);



%registration images
% %red 1040 aligned to calcium 920
% 
% reg_obj = imregcorr(red1040, calcium);
%add red1040+green 1040 
im1040=green+red1040;
% imALL=green+red1040+blue;
%% 
reg_obj = imregcorr(im1040, calcium);
% reg_obj = imregcorr(green, calcium);
T = reg_obj.T;

aligned_red1040 = imwarp(red1040, affine2d(T), 'OutputView', imref2d(size(red1040)));
imwrite(aligned_red1040, ([path,'aligned_red.tif']), 'tif');

aligned_green = imwarp(green, affine2d(T), 'OutputView', imref2d(size(green)));
imwrite(aligned_green, ([path,'aligned_green.tif']), 'tif');

% %red 890 aligned to red1040 and used to align Blue890 
reg_obj = imregcorr(red890, aligned_red1040);
T = reg_obj.T;
% aligned_red890 = imwarp(red890, affine2d(T), 'OutputView', imref2d(size(red890)));
% imwrite(aligned_red890, ([path,'aligned_red890.tif']), 'tif');

aligned_blue = imwarp(blue, affine2d(T), 'OutputView', imref2d(size(blue)));
imwrite(aligned_blue, ([path,'aligned_blue.tif']), 'tif');

%create composite image
composite(:,:,1)=aligned_red1040;
composite(:,:,2)=aligned_green;
composite(:,:,3)=aligned_blue;
%composite(:,:,4)=calcium;
imwrite(composite, ([path,'composite.tif']), 'tif');

imwrite(aligned_red1040, ([path,'compositecalcium.tif']), 'tif', 'WriteMode', 'overwrite');
imwrite(aligned_green, ([path,'compositecalcium.tif']), 'tif', 'WriteMode', 'append');
imwrite(aligned_blue, ([path,'compositecalcium.tif']), 'tif', 'WriteMode', 'append');
imwrite(calcium, ([path,'compositecalcium.tif']), 'tif', 'WriteMode', 'append');

disp ('done')