% nuc_segmentor_lite, will segment image with minimum algorithm for speed.
%
% [seg_img] = nuc_segmentor_lite3(img, nuc_rad, ss_factor, th_level, edge_option, refine_option)
%
% input: img: image to segment
%        nuc_rad: radius of nucleus, used to eliminate smaller object and
%        help object separation
%        ss_factor: X,Y subsampling. The higher the faster but the rougher
%        contours of the nucleus
%        opt_medif: median filter is used to get cleaner image before
%        segmenting. This is expensive in 3D, so option to turn it on. By
%        default, this filter is off. Recommended in 2D.
%        edge_option: 1 (default). IF positive, will removed nuclei
%        touching edges of the image
%        refile_option: 1 (default). Will refine each individual foci based
%        on initial segmentation (local threshold, only makes sense if
%        th_level not defined). This may drop nuclei that are close by.
%
% Separation of nuc based on distance transform.
% October 2009, Sylvain Costes, Lawrence Berkeley Lab
%

function [seg_img] = nuc_segmentor_lite3(img, nuc_rad, ss_factor, th_level,edge_option)

% Get image size
img = squeeze(img);
[w,h,d] = size(img);
%img = reshape(img,[w,h,d]); % make sure image is in 3D format for segmentation
dim = length(size(squeeze(img)));

%Subsample image only in X,Y plane
radS = round(nuc_rad/ss_factor);
if dim>2
    imgS = img(0:ss_factor:end,0:ss_factor:end,:);
else
    imgS = img(0:ss_factor:end,0:ss_factor:end);
end

% Use basic isodata threshold from Dipimage
if exist('th_level','var')
    if isempty(th_level)
        maskS = threshold(imgS);
    elseif isnan(th_level)
        maskS = threshold(imgS);
    else
        maskS = imgS > th_level;
    end
else
     maskS = threshold(imgS);
end

% Set edge_option to 1 if not entered. Default
if ~exist('edge_option','var')
    edge_option =1;
end

if edge_option
    if dim>2 % Get rid of any object on edges
        maskS(:,0,:) = 1; % Will get rid of anything on edge
        maskS(0,:,:) = 1;
        maskS(:,end,:) = 1;
        maskS(end,:,:) = 1;
    else
        maskS(:,0) = 1; % Will get rid of anything on edge
        maskS(0,:) = 1;
        maskS(:,end) = 1;
        maskS(end,:) = 1;
    end
    maskS = label(maskS>0) > 1;
end

% This will fill holes in mask
if dim>2
    for i=0:d-1
        inv_mask = label(~maskS(:,:,i));
        ms = measure(inv_mask,[],'size');
        [max_size,max_ind] = max(ms.size);
        max_ID = ms.ID(max_ind);
        maskS(:,:,i) = ~(inv_mask == max_ID); % Fill holes. Done for each slice, it works better.
    end
else
        inv_mask = label(~maskS);
        ms = measure(inv_mask,[],'size');
        [max_size,max_ind] = max(ms.size);
        max_ID = ms.ID(max_ind);
        maskS = ~(inv_mask == max_ID); % Fill holes. Done for each slice, it works better.
end
maskS = label(maskS,dim,round((pi*radS^dim)/2),100000000000000); % Remove small objects to clean up before spliting remaining
num_nuc = max(maskS);

% Use distance transform to separate touching nuclei
fprintf('Split object\n');
maskS = object_separation_2D(maskS>0,radS,imgS);

% Extrapolate back nuclear image
if ss_factor>1
    fprintf('Interpolating mask image back to full size\n');
    if dim>2
        seg_img = squeeze(oversample(maskS,[w,h,d]));
    else
        seg_img = squeeze(oversample(maskS,[w,h]));
    end
else
    seg_img = maskS;
end
seg_img = dip_image(uint16(seg_img));