function [label_mask] = object_separation_2D(mask,rad,gray_img)
% Same as object_separation, but will force a 3D image to be split solely
% based on 2D projection, then extrapolated back to 3D (usefull for cell
% monolayer)
%Routine that will separate ojbect from a binary mask based on their shape
%using distance transform
%   [label_mask] = object_separation(mask,rad,gray_img)
%   input:  mask, binary image defining objects
%           rad, radius of object
%           gray_img, optional. If in, will be used to weight seed growth
%   output: label_mask, binary mask with object separated and labeled
%
% Sylvain Costes, September 2009, Lawrence Berkeley Laboratory

if (max(mask) > 0) % make sure image is not blank
    mask = squeeze(mask>0); % make sure 2D image is not entered as 3D non binary
    dim = length(size(mask));
    if dim>2
        mask2D = label(squeeze(threshold(sum(mask,[],3))),2,10,1000000000)>0; % Remove small points that may appear with projection
    else
        mask2D = mask;
    end
    d = dt(mask2D); % Distance transform of binary mask
    min_mask = dip_localminima(-d,mask2D,2,100,rad,1); %size of local minima is factor that changes drastically the ouput.
    seed_img = bdilation(min_mask,round(rad/2));
    ms = measure(seed_img,[],'center');
    seed_img = newim(size(seed_img));
    seed_img(sub2ind([size(mask2D,2),size(mask2D,1)],round(ms.Center(2,:)),round(ms.Center(1,:)))-1) = ms.ID;
    if dim>2 % Need to grow back 2D seed image into a 3D image. Just place seeds on all slices where mask is present
        seed_img = repmat(seed_img,1,1,size(mask,3))*mask;
        d = repmat(d,1,1,size(mask,3));
    end
    seed_img = dip_image(uint16(seed_img));    
    if exist('gray_img')
        label_mask = dip_growregions(seed_img,d*gray_img,mask, dim, 0, 'high_first');
    else
        label_mask = dip_growregions(seed_img,[],mask, dim, 0, 'high_first');
    end
    % Remove nuclei smaller than half the desired size, and renumber nuclei with consecutive values
    % (tricky...)
    ms = measure(label_mask,[],'size');
    num_nuc = size(ms,1);
    temp = msr2obj(label_mask,ms,'size');
    temp= (temp>rad^dim*pi()/2)*mask;
    temp = dip_image(uint16(temp*label_mask));
    ms2 = measure(temp,[],'size');
    num_nuc = size(ms2,1);
    ms_new=dip_measurement(ms2.ID,'newID',1:num_nuc);
    label_mask = dip_image(uint16(msr2obj(temp,ms_new,'newID')));
else
    label_mask = mask;
end

