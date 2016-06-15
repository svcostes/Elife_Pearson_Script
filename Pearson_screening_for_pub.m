% This script looks into the target_dir and assumes all subdirectories have
% images and must be analyzed. Compute the
% correlation between the nuclear channel and two other channels within
% each nucleus. Nuclear segmentation is automatic. nuc_rad is the expected
% radius of the nucleus, if set too low, it will results in fragmented
% nuclei in the segmentation.
% S. Costes, January 2011, LBNL
% Version 2: compute 4 rings of equal bin, regardless of the nuclear size
%
%THIS SCRIPT WAS USED FOR THE GENOME-WIDE RNAi SCREEN (Swenson, J. et al., 2016, In prep)

%dipsetpref('NumberOfThreads',1); % Comment out if using 2D images
%Main directory
target_dir = 'example_images/';
output_dir = 'example_images/';

%%%%%%%%%%%%%%%%%%%%
file_ext = 'tif'; % THIS VESION ONLY READS TIF - CHANGE FOR YOUR ONW FILE FORMAT

%%%%%%%%% MANUAL SEGMENTING MODE: 1 = true, 0 = false
manual_threshold = 0;
font_size = 8;
line_thickness = 1; %minimum is 1, below, big PROBLEMS
background_correction = 1; % Set to 1 if want to correct nuclei background. Set to 0 if you dont want it
blurr_rad = 25; % the larger, the milder the background correction: 50 for the screen - Only usefull if previous option set to 1
disp_size = 100; % Size scale used to display image. 100% is the original size. If image is too large, use smaller value like 50 for 50%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
search_str = '*';
list_dir = dir([target_dir search_str]);
nuc_chan = 3; %for screen nuc_chan =3 (blue)
chan1 = 1; % Red fluorophore (channel 1 in rgb image)
chan2 = 2;  % Green fluorophore (channel 2 rgb image)
th_val = []; % IF set to [], then system does isodata to segment nucleus. If not, then uses constant threshold
nuc_rad = 15; % is the expected radius of the nucleus, if set too low, it will results in fragmented nuclei in the segmentation.

cnt=1;
for i=1:length(list_dir)
    if and(list_dir(i).isdir,~strcmp(list_dir(i).name(1),'.'))
        list_name{cnt} = list_dir(i).name;
        cnt = cnt+1;
    end
end
list_name = list_name(end:-1:1);

num_t = length(list_name);

for i_t=1:num_t % Open directory from list_name number i_t
    list_img = dir([target_dir list_name{i_t} '/*.' file_ext]);
    num_img = length(list_img);
    for i_img = 1:num_img % Read image i_img for time point i_t in directory list_name{i_t}
        fprintf('Trying to read %s\n',[target_dir list_name{i_t} '/' list_img(i_img).name]);
        temp = readim([target_dir list_name{i_t} '/' list_img(i_img).name],'TIFF'); % Read color tiff. Channel order is given above 
        name = list_img(i_img).name;
        if i_img == 1
            file_namej = ['Pearson_summary_',date]; %plate_name, date, pearson_summary_campus_40X
            ofp = fopen([output_dir file_namej '.txt'],'w');
            fprintf(ofp,'File name\tnuc#\tnuc size\tPearson with chan1\tPearson with chan2\tPearson chan1 & chan2\t');
            fprintf(ofp,'Shape factor\tsize ring1\tsize ring2\tsize ring3\tsize ring4\t');
            fprintf(ofp,'chan1 mean intensity\tchan1 std intensity\tchan1 Min intensity\tchan1 Max intensity\tchan1 Skew intensity\tchan1 Kurtosis intensity\t');
            fprintf(ofp,'chan1 mean ring1\tchan1 mean ring2\tchan1 mean ring3\tchan1 mean ring4\t');
            fprintf(ofp,'chan2 mean intensity\tchan2 std intensity\tchan2 Min intensity\tchan2 Max intensity\tchan2 Skew intensity\tchan2 Kurtosis intensity\t');
            fprintf(ofp,'chan2 mean ring1\tchan2 mean ring2\tchan2 mean ring3\tchan2 mean ring4\t');
            fprintf(ofp,'nuc mean intensity\tnuc std intensity\tnuc Min intensity\tnuc Max intensity\tnuc Skew intensity\tnuc Kurtosis intensity\t');
            fprintf(ofp,'nuc mean ring1\tnuc mean ring2\tnuc mean ring3\tnuc mean ring4\t');
        end
        try % some images are not with the same # of channels... Escape and print out there was a Pb with it
            % Add background correction
            nuc = temp{nuc_chan};
            if background_correction>0
                nuc = nuc - gaussf(nuc,blurr_rad);
            end
            if manual_threshold
                uiwait(dipshow(max(nuc,[],3),'lin'));
                th_val = str2double(inputdlg('Please input the manual defined threshold'));
            end
            % Do nuclear segmentation and remove deformed objects
            nuc_mask = nuc_segmentor_lite3(nuc,nuc_rad,2,th_val,1);
            ms = measure(nuc_mask,[],'p2a');
            nuc_mask = msr2obj(nuc_mask,ms,'p2a');
            nuc_mask = label(and(nuc_mask>0,nuc_mask<8)); % any object with a p2a greater than 8 are just too deformed to be a cell
            % Put blank pixel between touching nuclei (need this to have
            % ring computation not messed up between nuclei. 7/16/2012, SVC
            nuc_edge = dilation(nuc_mask,2)-nuc_mask;
            nuc_mask(nuc_edge>0) = 0;
            % create about 4 concentric rings around each nucleus
            nuc_ring = dt(nuc_mask>0); 
            ms_ring = measure(nuc_mask,nuc_ring,'MaxVal');
            ring_max = msr2obj(nuc_mask,ms_ring,'MaxVal');
            nuc_ring = round(nuc_ring*4./ring_max);
            nucr1 = nuc_ring==1;
            nucr2 = nuc_ring==2;
            nucr3 = nuc_ring==3;
            nucr4 = nuc_ring==4;
            num_nuc = max(nuc_mask);
 
            % add for mean intensity output
            ms1 = measure(nuc_mask,temp{chan1},{'mean','StdDev','Skewness','MaxVal','MinVal','ExcessKurtosis'});
            ms2 = measure(nuc_mask,temp{chan2},{'mean','StdDev','Skewness','MaxVal','MinVal','ExcessKurtosis'});
            msNuc = measure(nuc_mask, temp{nuc_chan},{'mean','StdDev','Skewness','MaxVal','MinVal','ExcessKurtosis','size','P2A'});
            ms1_mean = ms1.Mean; ms1_std = ms1.StdDev; ms1_sk = ms1.Skewness; ms1_min = ms1.MinVal; ms1_max = ms1.MaxVal; ms1_ex = ms1.ExcessKurtosis;
            ms2_mean = ms2.Mean; ms2_std = ms2.StdDev; ms2_sk = ms2.Skewness; ms2_min = ms2.MinVal; ms2_max = ms2.MaxVal; ms2_ex = ms2.ExcessKurtosis;
            n_mean = msNuc.Mean; n_size = msNuc.size; n_shape = msNuc.p2a; n_std = msNuc.StdDev; n_sk = msNuc.Skewness; n_min = msNuc.MinVal; n_max = msNuc.MaxVal; n_ex = msNuc.ExcessKurtosis;
            % measure concentric rings properties
            msr1 = measure(nuc_mask,nucr1*1.0,'sum');sizer1 = msr1.sum;
            msr2 = measure(nuc_mask,nucr2*1.0,'sum');sizer2 = msr2.sum;
            msr3 = measure(nuc_mask,nucr3*1.0,'sum');sizer3 = msr3.sum;
            msr4 = measure(nuc_mask,nucr4*1.0,'sum');sizer4 = msr4.sum;
               % measurement of Channel 1
            msr1 = measure(nuc_mask,nucr1*temp{chan1},'sum');sumr1_1=msr1.sum;meanr1_1=sumr1_1./sizer1;
            msr2 = measure(nuc_mask,nucr2*temp{chan1},'sum');sumr2_1=msr2.sum;meanr2_1=sumr2_1./sizer2;
            msr3 = measure(nuc_mask,nucr3*temp{chan1},'sum');sumr3_1=msr3.sum;meanr3_1=sumr3_1./sizer3;
            msr4 = measure(nuc_mask,nucr4*temp{chan1},'sum');sumr4_1=msr4.sum;meanr4_1=sumr4_1./sizer4;
               % measurement of Channel 2
            msr1 = measure(nuc_mask,nucr1*temp{chan2},'sum');sumr1_2=msr1.sum;meanr1_2=sumr1_2./sizer1;
            msr2 = measure(nuc_mask,nucr2*temp{chan2},'sum');sumr2_2=msr2.sum;meanr2_2=sumr2_2./sizer2;
            msr3 = measure(nuc_mask,nucr3*temp{chan2},'sum');sumr3_2=msr3.sum;meanr3_2=sumr3_2./sizer3;
            msr4 = measure(nuc_mask,nucr4*temp{chan2},'sum');sumr4_2=msr4.sum;meanr4_2=sumr4_2./sizer4;
               % measurement of nuclear channel
            msr1 = measure(nuc_mask,nucr1*temp{nuc_chan},'sum');sumr1_n=msr1.sum;meanr1_n=sumr1_n./sizer1;
            msr2 = measure(nuc_mask,nucr2*temp{nuc_chan},'sum');sumr2_n=msr2.sum;meanr2_n=sumr2_n./sizer2;
            msr3 = measure(nuc_mask,nucr3*temp{nuc_chan},'sum');sumr3_n=msr3.sum;meanr3_n=sumr3_n./sizer3;
            msr4 = measure(nuc_mask,nucr4*temp{nuc_chan},'sum');sumr4_n=msr4.sum;meanr4_n=sumr4_n./sizer4;
            % Measure pearson on per nucleus basis
            p1 =label_pearson_coef(temp{nuc_chan},temp{chan1},nuc_mask);
            p2 =label_pearson_coef(temp{nuc_chan},temp{chan2},nuc_mask);
            p1_2 =label_pearson_coef(temp{chan1},temp{chan2},nuc_mask);

            for i_nuc = 1:num_nuc
              fprintf(ofp,'%s\t%d\t%d\t%5.3f\t%5.3f\t%5.3f\t',name,i_nuc,n_size(i_nuc),p1(i_nuc),p2(i_nuc),p1_2(i_nuc));
              fprintf(ofp,'%5.3f\t%d\t%d\t%d\t%d\t',n_shape(i_nuc),sizer1(i_nuc),sizer2(i_nuc),sizer3(i_nuc),sizer4(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%d\t%d\t%5.3f\t%5.3f\t',ms1_mean(i_nuc),ms1_std(i_nuc),ms1_min(i_nuc),ms1_max(i_nuc),ms1_sk(i_nuc),ms1_ex(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%5.3f\t%5.3f\t',meanr1_1(i_nuc),meanr2_1(i_nuc),meanr3_1(i_nuc),meanr4_1(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%d\t%d\t%5.3f\t%5.3f\t',ms2_mean(i_nuc),ms2_std(i_nuc),ms2_min(i_nuc),ms2_max(i_nuc),ms2_sk(i_nuc),ms2_ex(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%5.3f\t%5.3f\t',meanr1_2(i_nuc),meanr2_2(i_nuc),meanr3_2(i_nuc),meanr4_2(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%d\t%d\t%5.3f\t%5.3f\t',n_mean(i_nuc),n_std(i_nuc),n_min(i_nuc),n_max(i_nuc),n_sk(i_nuc),n_ex(i_nuc));
              fprintf(ofp,'%5.3f\t%5.3f\t%5.3f\t%5.3f\t',meanr1_n(i_nuc),meanr2_n(i_nuc),meanr3_n(i_nuc),meanr4_n(i_nuc));
            end
            % Save nucleus segmentation for checking later
            % check if it's 3D
            if (size(nuc_mask,3) > 1)
                mask2D = squeeze(max(nuc_mask,[],3));
                nuc2D = squeeze(max(temp{nuc_chan},[],3));
            else
                mask2D = nuc_mask;
                nuc2D = temp{nuc_chan};
            end
            % Look up coordinates of nuc
            msn = measure(mask2D,[],{'Minimum'});
            mask2D = mask2D>0;
            overlay(stretch(nuc2D),mask2D-berosion(mask2D,line_thickness))
            % Label nuclei with their numbers
            for i_label=1:size(msn,1)
                text(msn(i_label).Minimum(1),msn(i_label).Minimum(2),sprintf('%d',i_label),'color','yellow','FontSize',font_size)
            end
            diptruesize(gcf,disp_size)
            saveas(gcf,[output_dir name(1:end-4) '_nuc_check.tif'],'tif');
            delete(gcf);
        catch
            fprintf(ofp,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',name,repmat(-1e6,8,1));
            fprintf('Problem with image #%d: %s\n',i_img,name);
        end
    end
    fclose('all');
end