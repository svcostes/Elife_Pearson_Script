function rp = label_pearson_coef(img1,img2,mask)
% label_pearson_coef(img1,img2,mask)
% Compute the pearson coefficient betweem img1 and img2 over mask region
% Each label in mask are treated separately. This is done without a loop
% using vectorization approaches and the measure function. So, it should
% compute for each nucleus in mask its pearson, fairly fast.
%
% Sylvain Costes, Berkeley Lab, May 2011

% Nomalize each nucleus
ms1 = measure(mask,img1,'mean');
ms2 = measure(mask,img2,'mean');
n1 = msr2obj(mask,ms1,'mean');
n2 = msr2obj(mask,ms2,'mean');
n1 = img1 - n1; % Mean subtracted nucleus
n2 = img2 - n2;
ms1 = measure(mask,n1^2,'sum');
ms2 = measure(mask,n2^2,'sum');
s1 = msr2obj(mask,ms1,'sum');
s2 = msr2obj(mask,ms2,'sum');
n1 = n1/sqrt(s1);
n2 = n2/sqrt(s2);

ms = measure(mask,n1*n2,'sum');
rp = ms.sum;