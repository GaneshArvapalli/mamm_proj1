 function [estdiag, estmaskleft,estmaskright] = runProject1(mammoimgleft,mammoimgright)
% Inputs:   mammoimgleft  -     the mammogram of the left side (rowL x colL)
%                               where rowL and colL are the image
%                               dimensions.
%           mammoimgright -     the mammogram of the right side (rowR x colR)
%
% Outputs:  estdiag      -   the estimated diagnosis.  Should only have
%                            values of 0(healthy), 1(benigh), and 2(cancer).
%                            and be of size (1 x 2) (left, right).
%
%           estmaskleft  -   the output binary mask for the left side. 
%                            Should only have values of zero or one.
%                            If the estdiag is 0, the mask should only have 
%                            values of 0. Should be of size (rowL x colL).
%           estmaskright -   the output binary mask for the right side. 
%                            Should only have values of zero or one.
%                            If the estdiag is 0, the mask should only have 
%                            values of 0. Should be of size (rowR x colR).
%


% Initializing variables
estdiag = zeros(1,2);
estmaskleft = zeros(size(mammoimgleft));
estmaskright = zeros(size(mammoimgright));

%% PUT IN YOUR DIAGNOSIS AND SEGMENTATION CODE BELOW!

% get intensity frequencies for both images
[l_counts,~] = imhist(mammoimgleft);
[r_counts,~] = imhist(mammoimgright);

% get sorted peak intensity frequencies for both image
[sort_lpk,sort_lpk_loc] = findpeaks(l_counts,'SortStr','descend');
[sort_rpk,sort_rpk_loc] = findpeaks(r_counts,'SortStr','descend');

% ignore any peak intensity frequencies that are black
sort_lpk(sort_lpk_loc < 50) = [];
sort_rpk(sort_rpk_loc < 50) = [];
sort_lpk_loc(sort_lpk_loc < 50) = [];
sort_rpk_loc(sort_rpk_loc < 50) = [];

% determine if the breast is healthy
if abs(sort_lpk(1)-sort_rpk(1))/min(sort_lpk(1),sort_rpk(1)) < 0.07
    return
end

% get peak intensity frequencies for both images
[lpk,lpk_loc] = findpeaks(l_counts);
[rpk,rpk_loc] = findpeaks(r_counts);

% ignore any peak intensity frequencies that are too dark
lpk_loc(lpk_loc < 100) = [];
rpk_loc(rpk_loc < 100) = [];

% determine which intensity the tumor is
[l_max,l_i] = max(lpk_loc([2:end])-lpk_loc([1:end-1]));
[r_max,r_i] = max(rpk_loc([2:end])-rpk_loc([1:end-1]));

% determine which breast has the cancer
lpk_m = mean(sort_lpk([1:5]));
rpk_m = mean(sort_rpk([1:5]));
left = 0;
if rpk_m > lpk_m
    left = 1;
end

% determine whether the tumor is malignant or a cancer
dis = abs(l_max-r_max)/min(l_max,r_max);
if dis > 0.7
    estdiag = [left,(1-left)];
else
    estdiag = [left*2,(1-left)*2];
end

% update the estimated mask
if left == 1
    estmaskleft(mammoimgleft > lpk_loc(l_i)) = 1;
else
    estmaskright(mammoimgright > rpk_loc(r_i)) = 1;
end
