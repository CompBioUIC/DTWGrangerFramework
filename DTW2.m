function [dist, warp] = DTW2(d1, d2, band_h, band_v)
% Developer: Ivan Brugere
% input:

%  d1 2-dimensional trajectory data of size [d x t] where d the x,y of data and t the time (e.g. [2x1000] matrix for trajectory of length 1000)
%  d2 2 dimensional trajectory data of size [d x t] where d the x,y of data and t the time
%  band_h the maximum fraction of warping to apply in the horizontal direction (e.g. 0.25)

%output:

%  dist the DTW distance between d1 and d2
%  warp the warping path through the DP matrix, represented as column j - row i

%% default band values
if (~exist('band_h','var') || isempty(band_h)),
    band_h=0.1;
end
if (~exist('band_v','var') || isempty(band_v)),
    band_v=0.1;
end

%size of matrix
[m,s] = size(d1);


%check d2 sizes
if(m ~= size(d2, 1) || s ~= size(d2, 2) || any(isnan(d1(:))) ||any(isnan(d2(:))))
    dist = NaN;
    warp = [];
    return;
elseif(m < s)
    d1 = d1';
    d2 = d2';
    [s,m] = size(d1);
end

    %% preallocate
    D = inf(s); % [s x s] DP matrix
    warp_path = NaN(s);  %[s x s] DP matrix for warping directionality
    
    %% band calculations and initial row column
    band_horizontal = band_h; 
    
    if (band_horizontal<=1)
        band_horizontal=ceil(band_h*s); 
    end;

    band_vertical = band_v;
    
    if (band_vertical<=1)
        band_vertical=ceil(band_v*s); 
    end;    
    
    band_lim_horizontal = min(1+band_horizontal,s);    
    d_temp = pdist2(d1(1, :), d2(1:band_lim_horizontal,: )); 
    D(1,1:band_lim_horizontal) = cumsum(d_temp);
    
    band_lim_vertical = min(1+band_vertical,s);    
    d_temp = pdist2(d2(1, :), d1(1:band_lim_vertical,: )); 
    D(1:band_lim_vertical,1) = cumsum(d_temp);
    
    
    %% main DP fill
    for i=2:s
        st = max(2,i-band_vertical);
        en = min(s,i+band_horizontal);
        %tic;
        for j= st:en
            [D(i,j), warp_path(i, j)] = min([D(i-1,j-1), D(i-1,j), D(i,j-1)]);            
            D(i,j)= D(i,j) + norm(d1(i, :) - d2(j, :));%+sum(sqrt(abs(d1(i, :) - d2(j, :)).^2)); 
        end
        %toc()
        
    end
    
    %% warping path reconstruction
    dist = D(s,s);
    i = s;
    j=s;
    k=s*2;
    warp = NaN(1, 2*s);
    while i ~= 1 && j ~= 1
        if(warp_path(i,j) == 1)
            i = i - 1;
            j = j - 1;
        elseif(warp_path(i,j) == 2)
            i = i - 1;
        elseif(warp_path(i,j) == 3)
            j = j - 1;
        end
        warp(k) = j - i;
        k = k - 1;
    end    
%     end
end
