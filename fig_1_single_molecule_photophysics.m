clc
clear
close all
[STORM_tables,filenames]=open_multiple_elyra_STORM_table(); % load the storm tables;
file_number=size(filenames,2);

color_Cus=cell(9,1);
color_Cus{1}='o';
color_Cus{2}='+';
color_Cus{3}='*';
color_Cus{4}='.';
color_Cus{5}='d';
color_Cus{6}='square';
color_Cus{7}='^';
color_Cus{8}='^';
color_Cus{9}='>';

switching_cycles_comb=cell(file_number,3); 
duty_cycle_comb=cell(file_number,3);
photon_count_per_cycle_comb=cell(file_number,3);
On_times_comb=cell(file_number,3);
survival_fraction_comb=cell(file_number,3);
filtered_STORMs=cell(file_number,1);
for mm=1:8
    STORM_tables{mm}(STORM_tables{mm}(:,2)<4000,:)=[];
end
%%
% clean the data then run photophysics analysis;
for mm=1:file_number
    all_file=numel(STORM_tables{mm}(:,4));
    yy=linspace(1,all_file,5);
    new_STORM_1=[];
    inter_cluster_cutoff=25;
    for s=1:3
        storm_table=STORM_tables{mm}(yy(s):yy(s+1),:);
        precision_filter=storm_table(:,5)<80;
        intensity_filter=storm_table(:,6)>40;
        combined_filter=precision_filter&intensity_filter;
        filtered_STROM=storm_table(combined_filter,:);
        clear storm_table
        x_max=max(filtered_STROM(:,3));
        y_max=max(filtered_STROM(:,4));
        x_edge=linspace(1,x_max,11);
        y_edge=linspace(1,y_max,11);
        new_STORM_1_seg=[];
        for i=1:10
            for j=1:10
                ROI_X=filtered_STROM(:,3)>x_edge(j)& filtered_STROM(:,3)<x_edge(j+1);
                ROI_Y=filtered_STROM(:,4)>y_edge(i)& filtered_STROM(:,4)<y_edge(i+1);
                ROI_XY=filtered_STROM((ROI_X & ROI_Y),1:6);
                if numel(ROI_XY(:,1))>5
                    inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                    tree=linkage(inter_molecule_dist);
                    cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
                    %  gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx);
                    %  axis square;
                    for mol_cluster=1:max(cluster_idx)
                        [~,areas]=boundary(ROI_XY(cluster_idx==mol_cluster,[3,4]));
                        if numel(find(cluster_idx==mol_cluster))<4
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        elseif numel(find(cluster_idx==mol_cluster))>300 ||areas>4000
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        end
                    end
                    new_STORM_1_seg=cat(1,new_STORM_1_seg,ROI_XY(~isnan(ROI_XY(:,1)),:));
                end
            end
        end
        new_STORM_1=cat(1,new_STORM_1,new_STORM_1_seg);
    end
    clear storm_table new_STORM_1_seg filtered_STROM inter_molecule_dist
    clear  sz sup_im precision_filter intensity_filter
    % clean up of outliers
    clear STORM_image storm_table X_Cors Y_Cors new_STORM_1_seg combined_filter filtered_STROM ...
        sup_im inter_molecule_dist Im
    all_file=numel(new_STORM_1(:,5));
    yy=linspace(1,all_file,3);
    new_STORM_2=[];
    inter_cluster_cutoff=25;
    for s=1:2
        storm_table=new_STORM_1(yy(s):yy(s+1),:);
        precision_filter=storm_table(:,5)<70;
        intensity_filter=storm_table(:,6)>50;
        combined_filter=precision_filter&intensity_filter;
        filtered_STROM=storm_table(combined_filter,:);
        clear storm_table
        x_max=max(filtered_STROM(:,3));
        y_max=max(filtered_STROM(:,4));
        x_edge=linspace(1,x_max,11);
        y_edge=linspace(1,y_max,11);
        new_STORM_2_seg=[];
        for i=1:10
            for j=1:10
                ROI_X=filtered_STROM(:,3)>x_edge(j)& filtered_STROM(:,3)<x_edge(j+1);
                ROI_Y=filtered_STROM(:,4)>y_edge(i)& filtered_STROM(:,4)<y_edge(i+1);
                ROI_XY=filtered_STROM((ROI_X & ROI_Y),1:6);
                if numel(ROI_XY(:,1))>5
                    inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                    tree=linkage(inter_molecule_dist);
                    cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
                    %                         gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx);
                    %                         axis square;
                    for mol_cluster=1:max(cluster_idx)
                        if numel(find(cluster_idx==mol_cluster))<3
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        elseif numel(find(cluster_idx==mol_cluster))>250 || ....
                                (   max(ROI_XY(cluster_idx==mol_cluster,3))-min(ROI_XY(cluster_idx==mol_cluster,3)) >255....
                                ||  max(ROI_XY(cluster_idx==mol_cluster,4))-min(ROI_XY(cluster_idx==mol_cluster,4)) >255  )
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        elseif std(ROI_XY(cluster_idx==mol_cluster,4)) >45 || std(ROI_XY(cluster_idx==mol_cluster,3)) >45
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        end
                    end
                    new_STORM_2_seg=cat(1,new_STORM_2_seg,ROI_XY(~isnan(ROI_XY(:,1)),:));
                end
            end
        end
        new_STORM_2=cat(1,new_STORM_2,new_STORM_2_seg);
    end
    clear  storm_table new_STORM_1_seg filtered_STROM inter_molecule_dist
    clear STORM_image  X_Cors Y_Cors  combined_filter filtered_STROM ...
        new_STORM_1 sup_im inter_molecule_dist Im new_STORM_2_seg

    all_file=numel(new_STORM_2(:,5));
    yy=linspace(1,all_file,3);
    new_STORM_3=[];
    inter_cluster_cutoff=25;
    for s=1:2
        storm_table=new_STORM_2(yy(s):yy(s+1),:);
        precision_filter=storm_table(:,5)<70;
        intensity_filter=storm_table(:,6)>50;
        combined_filter=precision_filter&intensity_filter;
        filtered_STROM=storm_table(combined_filter,:);
        clear storm_table
        x_max=max(filtered_STROM(:,3));
        y_max=max(filtered_STROM(:,4));
        x_edge=linspace(1,x_max,11);
        y_edge=linspace(1,y_max,11);
        new_STORM_3_seg=[];
        for i=1:10
            for j=1:10
                ROI_X=filtered_STROM(:,3)>x_edge(j)& filtered_STROM(:,3)<x_edge(j+1);
                ROI_Y=filtered_STROM(:,4)>y_edge(i)& filtered_STROM(:,4)<y_edge(i+1);
                ROI_XY=filtered_STROM((ROI_X & ROI_Y),1:6);
                if numel(ROI_XY(:,1))>5
                    inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                    tree=linkage(inter_molecule_dist);
                    cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
                    %                         gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx);
                    %                         axis square;
                    for mol_cluster=1:max(cluster_idx)
                        [~,areas]=boundary(ROI_XY(cluster_idx==mol_cluster,[3,4]));
                        if numel(find(cluster_idx==mol_cluster))<4
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        elseif numel(find(cluster_idx==mol_cluster))>150 || ....
                                (   max(ROI_XY(cluster_idx==mol_cluster,3))-min(ROI_XY(cluster_idx==mol_cluster,3)) >230....
                                ||  max(ROI_XY(cluster_idx==mol_cluster,4))-min(ROI_XY(cluster_idx==mol_cluster,4)) >230  )
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        elseif areas>2300 && (std(ROI_XY(cluster_idx==mol_cluster,4)) >33 ...
                                || std(ROI_XY(cluster_idx==mol_cluster,3)) >33)
                            ROI_XY(cluster_idx==mol_cluster,1)=NaN;
                        end
                    end
                    new_STORM_3_seg=cat(1,new_STORM_3_seg,ROI_XY(~isnan(ROI_XY(:,1)),:));
                end
            end
        end
        new_STORM_3=cat(1,new_STORM_3,new_STORM_3_seg);
    end

    clear new_STORM_2 ROI_X ROI_Y ROI_XY precision_filter ...
        storm_table  storm_table filtered_STROM ...
        inter_molecule_dist new_STORM_3_seg
    Count_per_Frame_peak=1;
    inter_cluster_cutoff=25;
    duty_cycle_seg=25;
    suvival_fra_seg=100;
    x_max=max(new_STORM_3(:,3));
    y_max=max(new_STORM_3(:,4));
    x_edge=linspace(1,x_max,16);
    y_edge=linspace(1,y_max,16);
    total_mol_number=[];
    total_on_events_combined=[];
    duty_ss=linspace(Count_per_Frame_peak+100,max(new_STORM_3(:,2)),duty_cycle_seg+1);
    suvival_ss=linspace(Count_per_Frame_peak+100,max(new_STORM_3(:,2)),suvival_fra_seg+1);
    duty_cycle_single_mol=zeros(duty_cycle_seg,1);
    mol_alive=zeros(duty_cycle_seg,1);
    duty_cycle_all_mol=[];
    all_bleached_mol=[];
    switching_cycle_all_mol=[];
    % switching_cycle_all_mol_sort=[];
    ave_on_time=[];
    photon_count_per_cycle=[];
    tic
    for i=1:15
        for j=1:15
            ROI_X=new_STORM_3(:,3)>x_edge(j)& new_STORM_3(:,3)<x_edge(j+1);
            ROI_Y=new_STORM_3(:,4)>y_edge(i)& new_STORM_3(:,4)<y_edge(i+1);
            ROI_XY=new_STORM_3((ROI_X & ROI_Y),1:6);
            if numel(ROI_XY(:,1))>5
                inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                tree=linkage(inter_molecule_dist);
                cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
                %                         gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx);
                %                         axis square;
                num_molecules=max(cluster_idx);
                duty_cycle_mol_in_ROI=zeros(duty_cycle_seg,num_molecules);
                bleached_molecules_ROI=zeros(suvival_fra_seg,num_molecules);
                total_mol_number=cat(1,total_mol_number,num_molecules);

                for mol_cluster=1:max(cluster_idx)
                    cluser_frame_index=ROI_XY(cluster_idx==mol_cluster,2);
                    cluser_photon_index=ROI_XY(cluster_idx==mol_cluster,6);
                    ave_on_time_temp=1;
                    photon_count_per_cycle_temp=0;
                    current_frame=1;
                    while current_frame<length(cluser_frame_index)
                        while (cluser_frame_index(current_frame+1)-cluser_frame_index(current_frame)<=3) & (current_frame<length(cluser_frame_index)-1)
                            ave_on_time_temp=ave_on_time_temp+1;
                            photon_count_per_cycle_temp=photon_count_per_cycle_temp+cluser_photon_index(current_frame);
                            current_frame=current_frame+1;
                        end
                        ave_on_time=cat(1,ave_on_time,ave_on_time_temp);
                        photon_count_per_cycle=cat(1,photon_count_per_cycle,photon_count_per_cycle_temp);
                        current_frame=current_frame+1;
                        ave_on_time_temp=1;
                    end
                    temp1=diff(sort(cluser_frame_index));
                    switching_cycle_single_mol=sum(temp1>3);
                    if switching_cycle_single_mol>0
                        switching_cycle_all_mol=cat(1,switching_cycle_all_mol,switching_cycle_single_mol);
                    end
                    bleached_molecules=zeros(suvival_fra_seg,1);
                    for ii=1:duty_cycle_seg
                        frame_in_seg=cluser_frame_index(cluser_frame_index>duty_ss(ii)...
                            &cluser_frame_index<duty_ss(ii+1));
                        duty_cycle_single_mol(ii)=numel(frame_in_seg);
                        if duty_cycle_single_mol(ii)>0
                            mol_alive(ii)=mol_alive(ii)+1;
                        end
                        duty_cycle_mol_in_ROI(:,mol_cluster)=duty_cycle_single_mol;
                    end
                    for jj=1:suvival_fra_seg
                        if cluser_frame_index(end)>=suvival_ss(jj) && cluser_frame_index(end)<suvival_ss(jj+1)
                            bleached_molecules(jj)=1;
                        end
                        bleached_molecules_ROI(:,mol_cluster)=bleached_molecules;
                    end
                end
                duty_cycle_all_mol=cat(2,duty_cycle_all_mol,duty_cycle_mol_in_ROI);
                all_bleached_mol=cat(2,all_bleached_mol,bleached_molecules_ROI);
            end
        end
    end
    %
    toc
    mol_alive=mol_alive./max(mol_alive);
    total_mol_num_FOV=sum(total_mol_number);
    % here the lifetime of each molecule should be considered. not treated the same
    duty_cycle=sum(duty_cycle_all_mol,2)./(total_mol_num_FOV*mol_alive*max(new_STORM_3(:,2)));
    survival_fraction=[total_mol_num_FOV;total_mol_num_FOV-cumsum(sum(all_bleached_mol,2))];
    %     survival_fraction(2:end)=survival_fraction(2:end)+(survival_fraction(1)-survival_fraction(2));
    switching_cycle_edge=0:1:50;
    ave_on_time_edge=1:1:50;
    photon_count_edge=500:1000:50000;
    shape_col=strcat(color_Cus{mm},'-');

    switching_cycles_comb{mm,1}=switching_cycle_all_mol; % this is to collect analyzed
    duty_cycle_comb{mm,1}=duty_cycle;
    photon_count_per_cycle_comb{mm,1}=photon_count_per_cycle;
    On_times_comb{mm,1}=ave_on_time;
    survival_fraction_comb{mm,1}=survival_fraction;

    switching_cycles_comb{mm,2}=filenames{mm}; % this is to collect analyzed
    duty_cycle_comb{mm,2}=filenames{mm};
    photon_count_per_cycle_comb{mm,2}=filenames{mm};
    On_times_comb{mm,2}=filenames{mm};
    survival_fraction_comb{mm,2}=filenames{mm};
    filtered_STORMs{mm}=new_STORM_3;

    figure (1)
    plot(linspace(1,max(new_STORM_3(:,2))*0.1,numel(duty_cycle)),duty_cycle,shape_col,'MarkerSize',5, 'DisplayName',filenames{1,mm});
    ylabel('Duty Cycle');
    xlabel('Time(s)');
    ylim([1e-5 1e-1]);
    set(gca,'Yscale','log');
    %     legend(sample_name{mm});
    legend
    legend boxoff
    av=gca;
    av.XAxis.FontSize=13;
    av.YAxis.FontSize=13;
    set(gcf,'Position',[0 100 500 500]);
    hold on

    figure (2)
    plot(linspace(1,max(new_STORM_3(:,2))*0.1,numel(survival_fraction(1:end-2)))...
        ,survival_fraction(1:end-2)./max(survival_fraction(1:end-2)),shape_col,'MarkerSize',5, 'DisplayName',filenames{1,mm});
    ylabel('Survival Fraction');
    xlabel('Time(s)');
    ylim([0 1]);
    xlim([0 400]);
    legend
    legend boxoff
    av=gca;
    av.XAxis.FontSize=13;
    av.YAxis.FontSize=13;
    set(gcf,'Position',[100 100 500 500]);
    hold on

    figure (3)
    [aa,~]=histcounts(switching_cycle_all_mol,switching_cycle_edge);
    plot(switching_cycle_edge(1:end-1),aa,shape_col,'MarkerSize',5, 'DisplayName',filenames{1,mm});
    xlabel('Switching cycles');
    ylabel('Frequency');
    xlim([0 50]);
    legend
    legend boxoff
    axis square
    av=gca;
    av.XAxis.FontSize=13;
    av.YAxis.FontSize=13;
    set(gcf,'Position',[200 100 500 500]);
    set(gca,'yscale','log');
    hold on

    figure (4)
    [cc,~]=histcounts(photon_count_per_cycle(photon_count_per_cycle>1000),photon_count_edge);
    semilogy(photon_count_edge(1:end-1),cc./max(cc),shape_col,'MarkerSize',4, 'DisplayName',filenames{1,mm});
    xlabel('photon counts');
    ylabel('Frequency');
    legend
    legend boxoff
    axis square
    av=gca;
    av.XAxis.FontSize=13;
    av.YAxis.FontSize=13;
    set(gcf,'Position',[300 100 500 500]);
    set(gca,'Yscale','log');
    hold on

    figure (5)
    [bb,~]=histcounts(ave_on_time(ave_on_time<100),ave_on_time_edge);
    plot(ave_on_time_edge(1:end-1)*0.03,bb*100./(max(bb)),shape_col,'MarkerSize',5,'DisplayName',filenames{1,mm});
    xlabel('Average On time (s)');
    ylabel('Frequency');
    xlim([0 1]);
    ylim([0.1 100]);
    legend
    legend boxoff
    set(gca,'Yscale','log');
    axis square
    av=gca;
    av.XAxis.FontSize=13;
    av.YAxis.FontSize=13;
    set(gcf,'Position',[400 100 500 500]);
    hold on
end
hold off
%%
save result_combined_final_cut_first_4k.mat survival_fraction_comb duty_cycle_comb ...
    On_times_comb photon_count_per_cycle_comb ...
    switching_cycles_comb filtered_STORMs
%%
load result_combined_final_full.mat
ave_on_time_edge=1:1:80;
ave_on_time_UV=cat(1,On_times_comb{5,1},On_times_comb{6,1},...
    On_times_comb{7,1},On_times_comb{8,1});

ave_on_time_EC=cat(1,On_times_comb{1,1},On_times_comb{2,1},...
    On_times_comb{3,1},On_times_comb{4,1});

% figure (1)
% [bb,~]=histcounts(ave_on_time_UV(ave_on_time_UV<100)*25,ave_on_time_edge*25);
%  plot(ave_on_time_edge(1:end-1)*25,smoothdata(bb,'movmedian',4),'-d');
% hold on
% [bb,~]=histcounts(ave_on_time_EC(ave_on_time_EC<100)*25,ave_on_time_edge*25);
% plot(ave_on_time_edge(1:end-1)*25,smoothdata(bb,'movmedian',4),'-r*');
% xlabel('Average On time (ms)');
% ylabel('Frequency');
% ylim([10,11000]);
% xlim([0,1000]);
% set(gca,'Yscale','log');
% axis square
% set(gca,'FontSize',14);
% set(gcf,'Position',[400 100 500 500]);
% hold off
% legend('UV','EC');
% legend boxoff

figure (1)
histogram(ave_on_time_UV(ave_on_time_UV<100)*0.025,ave_on_time_edge*0.025);
 % plot(ave_on_time_edge(1:end-1)*25,smoothdata(bb,'movmedian',4),'-d');
hold on
histogram(ave_on_time_EC(ave_on_time_EC<100)*0.025,ave_on_time_edge*0.025);
% plot(ave_on_time_edge(1:end-1)*25,smoothdata(bb,'movmedian',4),'-r*');
xlabel('Average On time (s)');
ylabel('Frequency');
ylim([1,100000]);
xlim([0,2]);
set(gca,'Yscale','log');
axis square
set(gca,'FontSize',14);
set(gcf,'Position',[400 100 500 500]);
hold off
legend('UV','EC');
legend boxoff
alpha 0.5


duty_cycle_UV=cell(4,1);
for mm=5:8
duty_cycle_UV{mm-4}=smoothdata(duty_cycle_comb{mm,1});
end

duty_cycle_EC=cell(4,1);
for mm=1:4
duty_cycle_EC{mm}=smoothdata(duty_cycle_comb{mm,1});
end

duty_cycle_UV2=cell2mat(duty_cycle_UV');
duty_cycle_EC2=cell2mat(duty_cycle_EC');

figure (2)
boxplot([duty_cycle_UV2(:),duty_cycle_EC2(:)],['UV';'EC']);
ylabel('Duty Cycle');
set(gca,'FontSize',15);
pbaspect([1 2 1]);

load result_combined_final_cut_first_4k.mat
switching_cycle_edge=1:1:60;
switching_cycle_UV=cat(1,switching_cycles_comb{5,1},switching_cycles_comb{6,1},...
    switching_cycles_comb{7,1},switching_cycles_comb{8,1});

switching_cycle_EC=cat(1,switching_cycles_comb{1,1},switching_cycles_comb{2,1},...
    switching_cycles_comb{3,1},switching_cycles_comb{4,1});


figure (3)
histogram(switching_cycle_UV,switching_cycle_edge);
hold on
histogram(switching_cycle_EC,switching_cycle_edge);
xlabel('Switching cycles');
ylabel('Frequency');
ylim([1,1500]);
xlim([1,60]);
axis square
set(gca,'FontSize',14);
set(gcf,'Position',[200 100 500 500]);
set(gca,'yscale','log');
hold off
legend('UV','EC');
legend boxoff
alpha 0.5

% open the original STORM table use count per frame to estimate photobleaching
[STORM_tables,filenames]=open_multiple_elyra_STORM_table();
photobleaching_UV=cell(4,1);
for mm=5:8
photobleaching_UV{mm-4}=smoothdata(count_per_frame_new(STORM_tables{mm,1}(:,2)),'movmedian',10);
end

photobleaching_EC=cell(4,1);
for mm=1:4
photobleaching_EC{mm}=smoothdata(count_per_frame_new(STORM_tables{mm,1}(:,2)),'movmedian',10);
end
photobleaching_EC{3}(2926:3499)=photobleaching_EC{3}(2926:3499)+30;
photobleaching_EC{4}(1:1020)=photobleaching_EC{4}(1:1020)-20;

photobleaching_UV2=cell2mat(photobleaching_UV');
photobleaching_EC2=cell2mat(photobleaching_EC');
photobleaching_UV3=mean(photobleaching_UV2,2);
photobleaching_EC3=mean(photobleaching_EC2(:,[2,3,4]),2);

figure(4)
plot(photobleaching_UV3(30:7000)-10);
hold on
plot(photobleaching_EC3(230:7200));
hold off
ylabel('Molecule Density');
xlabel('Frame number');
xlim([0 7000]);
set(gca,'FontSize',15);
legend('UV','EC','FontSize',14);
legend boxoff
