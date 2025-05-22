%% load the thunderstorm produced STORM table for moelcule density plot;
one=readmatrix('1 Hz.csv');
five=readmatrix('5 Hz.csv');
ten=readmatrix('10 hz.csv');
twenty=readmatrix('20 hz.csv');

one_density=count_per_frame(one(:,2));
five_density=count_per_frame(five(:,2));
ten_density=count_per_frame(ten(:,2));
twenty_density=count_per_frame(twenty(:,2));

x=0.02:0.02:0.02*5000;
plot(x-42,one_density-50,'LineWidth',1.2);
hold on
plot(x-42,five_density,'LineWidth',1.2);
plot(x-42,ten_density+50,'LineWidth',1.2);
plot(x-42,twenty_density+130,'LineWidth',1.2);
hold off
xlim([0 15]);
ylim([170 770]);
legend('1 Hz','5 Hz','10 Hz','20 Hz','NumColumns',4);
legend boxoff
xlabel('Time (s)');
ylabel('Molecule Density');
set(gca,'FontSize',14);

%% load the raw CZI file
clear
one=convert_CZI_to_tiff();
five=convert_CZI_to_tiff();
ten=convert_CZI_to_tiff();
twenty=convert_CZI_to_tiff();
combined={one,five,ten,twenty};
clear one five ten twenty 
%%
raw_image_mean=zeros(5000,4);
for i=1:4
    raw_image_mean(:,i)=squeeze(mean(combined{i},[1 2]));
end

x=(1:500)*0.022;
for i=1:4
plot(x,raw_image_mean(980:1479,i)./min(raw_image_mean(980:1479,i)),'LineWidth',1.5);
hold on
end
hold off
xlabel('Time(s)');
ylabel('Normalized Intensity');
xlim([0 10]);
set(gca,'FontSize',15);
legend('1 Hz','5 Hz','10 Hz','20 Hz','NumColumns',4);
legend boxoff

%% do pixel wised fft ananlysis.
fft_pix=zeros(400,404,2500);
ifft_pix=zeros(400,404,2500);
fft_result_combined=cell(4,1);
% aa=mode(cell2mat(combined),'all');

for ss=1:4
    for i=1:404
        for j=1:400
            pix=squeeze(combined{ss}(j,i,1:2500));
            pix=detrend(double(pix));
            pix_centered=pix-mean(pix);
            fft_pix(j,i,:)=fft(pix_centered);
        end
    end

    for i=1:404
        for j=1:400
            ifft_pix(j,i,:)=abs(fftshift(fft_pix(j,i,:)));
        end
    end
    fft_result_combined{ss}=ifft_pix;
end
clear fft_pix ifft_pix
%%
fft_result_combined_cut=fft_result_combined;
clear fft_result_combined
aa=127; % threshold to remove background, and set them as NaN;
for ss=1:4
    bb=mean(combined{ss}(1:400,1:404,:),3);
    cc=bb<aa;
    for i=1:2500
    dd=fft_result_combined_cut{ss}(:,:,i);
    dd(cc)=NaN;
    fft_result_combined_cut{ss}(:,:,i)=dd;
    end
end


t=22.1/1000; % t is the frame time;
fs=1/t; % aquisition frequency;
f=linspace(-fs/2,fs/2,2500); 
fft_1hz_max=mean(fft_result_combined_cut{1},[1,2],'omitnan');
fft_1hz_max=squeeze(fft_1hz_max);
fft_1hz_max=max(fft_1hz_max(f>0.5));

figure (1)
for ss=[1,2,3,4]
    spectrum=mean(fft_result_combined_cut{ss},[1,2],'omitnan');
    plot(f,squeeze(spectrum),'LineWidth',2);
    hold on
end
hold off
xlim([-20.5,20.5]);
xlabel('Frequency');
ylabel('Amplitude');
set(gca,'FontSize',14);
legend('1 Hz','5 Hz','10 Hz','20 Hz','FontSize',14);
legend boxoff


