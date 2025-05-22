% load the raw czi file and split it into x sized subsequences, calculate
% sofi for each subsequences, do drift correction then sum the images;
clear
close all
[files,root_filepath]=uigetfile('*.czi',MultiSelect='on');
file_number=length(files);
sub_seq_size=200;
for file_ind=1:file_number
    file=bfopen(files{file_ind});
    movie=file{1}(:,1);
    clear file
    frame_mean=zeros(size(movie,1),1);
    for frame_ind=1:size(movie,1)
        frame_mean(frame_ind)=mean(movie{frame_ind},[1 2]);
    end
    spike=find(frame_mean<0.65*mode(frame_mean)|...
        frame_mean>1.3*mode(frame_mean));
    movie(spike)=[];
    sub_folder_name=sprintf(files{file_ind}(1:end-4));
    new_folder_fullpath=strcat(root_filepath,sub_folder_name);
    mkdir(new_folder_fullpath);
    cd(new_folder_fullpath);
    stack=1;
    j=1;
    while j<=length(movie)
        if rem(j,sub_seq_size)==0
            stack=stack+1;
        end
        name=sprintf('stack_%d.tiff',stack);
        imwrite(mat2gray(movie{j}),name,"WriteMode","append");
        j=j+1;
    end
    clear movie
    temp=dir('*.tiff');
    stack_number=length(temp);

    for i=1:stack_number
        sofi_raw_stack{i}=temp(i).name;
        sofi_tif_size=numel(imfinfo(sofi_raw_stack{i}));
        if sofi_tif_size<50
            delete (sofi_raw_stack{i}); % remove the small sofi_tiff stacks
        end
    end
    temp=dir('*.tiff');
    stack_number=length(temp);
    img={stack_number};
    sofi_raw_stack={stack_number};
    sofi_stack=cell(length(temp),1);

    for i=1:stack_number
        sofi_raw_stack{i}=temp(i).name;
    end

    for i=1:stack_number
        sofi_tif_size=numel(imfinfo(sofi_raw_stack{i}));
        [sofi,grid]=sofiCumulants(sofi_raw_stack{i},1,sofi_tif_size);
        [sofi,fwhm]=sofiFlatten([],sofi,grid);
        sofi=sofiLinearize(sofi,fwhm);% "blind" linearization (no parameters)
        sofi_stack{i}=sofi;
        delete(sofi_raw_stack{i});
    end

    % make a new folder to store the calculated sofi image
    currentFolder = pwd;
    new_file_path=strcat(currentFolder,'\generated_sofi_image');
    mkdir(new_file_path);
    cd(new_file_path);

    % perform drift correction for each sofi order.
    for order=2:6
        order_name=sprintf('order_%d',order);
        sofi_folder_path=strcat(new_file_path,'\',order_name);
        mkdir(sofi_folder_path);
        cd(sofi_folder_path);
        ref_im=sofi_stack{1,1}{1,order};
        fft_one=fft2(ref_im);
        sofi_sum=ref_im;
        for i=2:numel(sofi_stack)
            target_im=sofi_stack{i,1}{1,order};
            fft_two=fft2(target_im);
            [shift, ~] = dftregistration(fft_one,fft_two,20);
            drift_corrected_im=imtranslate(target_im,[shift(4),shift(3)]);
            sofi_sum=sofi_sum+drift_corrected_im;
            stack_name=sprintf('_drift_cor_%d.tif',i);
            drift_corrected_im_name=strcat(order_name,stack_name);
            imwrite(mat2gray(drift_corrected_im),drift_corrected_im_name);
        end
        first_frame_im_name=strcat(order_name,'_ref_1.tif');
        sofi_sum_name=strcat(order_name,'sofi_sum.tif');
        imwrite(mat2gray(ref_im),first_frame_im_name);
        cd(new_file_path);
        imwrite(mat2gray(sofi_sum),sofi_sum_name);
    end
    cd(root_filepath);
end

%% do frc analysis;
FSC_SOFI=cell(5,1);
frequency=cell(5,1);
for i=1:5
current_folder=pwd;
[file_names,filepath]=uigetfile('*.tif','MultiSelect','on');
cd (filepath);
image_number=numel(file_names);
SOFI_1=double(imread(file_names{1}));
im_high=numel(SOFI_1(:,1));
im_width=numel(SOFI_1(1,:));
SOFI_1=SOFI_1(1:im_high,1:im_high);
sofi_1st_half=zeros(im_high,im_width);
sofi_2nd_half=zeros(im_high,im_width);
mid_split=image_number/2;
for j=1:2:image_number-1
    sofi_1st_half=sofi_1st_half+double(imread(file_names{j}));
    sofi_2nd_half=sofi_2nd_half+double(imread(file_names{j+1}));
end
gussian_width=round(sqrt(numel(SOFI_1))*0.2);
dx=gussian_width;
dy=gussian_width;
gaussian_window_size=(sqrt(numel(SOFI_1))-1)/2;

[X,Y] = meshgrid(-gaussian_window_size:gaussian_window_size); % overall PSF size;
fK = exp(-X.^2/dx^2-Y.^2/dy^2); % simulate PSF;
fK = fK/sum(fK(:));
fK_nom=rescale(fK,0,1);
order=i+1;

[FSC_SOFI{i},frequency{i}]=FSC_SOFI_2024(double(sofi_1st_half).*fK_nom,double(sofi_2nd_half).*fK_nom,order);
cd (current_folder);
plot(frequency{i},FSC_SOFI{i});
hold on 
end
hold off
legend('SOFI-2','SOFI-3','SOFI-4','SOFI-5','SOFI-6');

%% frc code
function  [FSC_final,frequency_range]=FSC_SOFI_2024(img1,img2,order)
dimension=length(img1(:,1));
pixel_size=100/order; 
max_freq=dimension;
im1_fft =fft2(img1,max_freq,max_freq);
im2_fft=fft2(img2,max_freq,max_freq);
cent_ima1_fft=abs(fftshift(im1_fft));
cent_ima2_fft=abs(fftshift(im2_fft));
[columnsInImage, rowsInImage] = meshgrid(1:max_freq, 1:max_freq);
centerX = max_freq/2;
centerY = max_freq/2;
ring_number=80;
FSC_final=zeros(ring_number,1);
highest_frequency=1/(pixel_size);
frequency_range=linspace(0,highest_frequency,ring_number);
ring_increment_step=ceil((max_freq/2)/ring_number);
for m=1:ring_number
    inner_circle = (-ring_increment_step+m*ring_increment_step)*1;
    outer_circle=ceil((m*ring_increment_step)*1);
    
    circlePixels = (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 <= outer_circle.^2 & (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 >= inner_circle.^2 ;
    lin_A=cent_ima1_fft(:);
    lin_B=cent_ima2_fft(:);
    ring=find(circlePixels==1 );
    delta_A=lin_A(ring)-mean(lin_A(ring));
    delta_B=lin_B(ring)-mean(lin_B(ring));
    FSC_final(m)=(delta_A'*delta_B)/sqrt(sum(delta_A.^2)*sum(delta_B.^2));
end
end
