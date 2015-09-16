function [D_sum, D_max] = dt_distance(samples_model, samples_exp, ...
    directions, stat1d, smoothing) 
% Wrapper for computing some "distances" between two datasets of the 
% Lignum model.

% INPUT:
% samples_model - size dim*n, samples from the model, each column == one
% datapoint
% samples_exp - size dim*m, experimental data corresponding the model, each 
% column == one datapoint
% directions - number of directions (projections) for which the 1d statistics 
% are computed (default 100)
% stat1d - some distance between the data sets (projection distance idea is
% used to define the 1d distances below to multiple dimensional case. 
% 1d distances are those with stat1d > 0)
% 1 == Kolmogorov-Smirnov i.e. sup-distance (default)
% 2 == Cramer von Mises i.e. L_2 distance
% 3 == Anderson-Darling, weighted L_2 distance
% 4 == Watson, modified Cramer von Mises distance
% 5 == L^1, L_1 distance
% 6 == L^2, basically the same as Cramer von Mises
% 7 == weighted L^2 basically the same as Anderson-Darling, 
% 8 == ecdf area^2, the area between empirical cdf's squared
% 9 == ecdf area, the area between empirical cdf's
% 10 == 1d energy distance
% -1 == NN distance
% -2 == Distance based on empirical moments of data
% -3 == Energy distance
% -4 == chi^2 type NN distance
% smoothing - 0 no smoothing, 1 if smooth the empirical cdfs (this makes the 
% distance continuous function which can useful for some purposes), 2 if 
% smooth using interpolation (default 0)
%
% OUTPUT:
% D_sum - average of 1d statistics computed to different directions (when stat1d > 0)
% D_max - max value of 1d statistics computed to different directions (when stat1d > 0)
%
% NOTE1:
% There are many other distances that could be used. The proper distance
% to choose depends on many settings i.e. what characteristics of the
% densities should be matched (e.g. are the tails important).
%
% NOTE2:
% Many other distances that are not implemented here could suit better; e.g.
% Kullback-Leibler and Hellinger are such (however they are rather difficult 
% to compute).
%
% NOTE3:
% These codes are for testing use and originally not intended for the 
% current use. Use with your own risk! Some distance methods could be also
% further tuned or could be made faster.
%
% NOTE4:
% About recommended settings: 
% 1) for projection pursuit distances, set directions to atleast 100; the
% higher the number the better D_sum and D_max approximate the projection pursuit 
% distance functions but, unfortunately, computational cost increases 
% linearly on the number of directions
% 2) Kolmogorov-Smirnov or Cramer von Mises are probably the best choices but 
% all distances should work; they just compare the datasets somewhat 
% differently; from the non-projection type distances the energy distance 
% is perhaps the most useful and it has nice theory behind it, unfortunately 
% it gets very slow if the dataset is large 
% 3) smoothing 2 seems to be  the best option for approximating the exact 
% value of the distance (smoothing 0 == exact value) with continuous function 
% when minimizing the distance with gradient-based methods. If one simply 
% wants to compute the distance then there is no reason to set smoothing to
% other value than 0. 
% 4) D_sum is the recommended option for the distance measure, d_max can be
% somewhat inaccurate unless the amount of directions is high, like >1000. 
%
% NOTE5:
% When comparing multiple data sets, then this function must be called
% several times. Some care must be taken how to weight each (sub)distance, 
% especially if the number of data varies. Some distances depend on scaling
% and some the amount of data.

% --- some settings/fixed values:
two_sample = 1; % 1 if compare two finite sets of data, 0 if model cdf is known
d1_scaling = 1; % additional scaling
d2_scaling = 1; % additional scaling
max_directions = 1000; % upper bound for the number of the directions
use_normalised_stat = 1; % 1 no scaling, 2 scaling ~ sqrt(n) or ~ sqrt(n1*n2/(n1+n2)) 
% so that the value of the distance does not approximately depend on sample size
% ---

% check input
dim = size(samples_model,1);
if nargin < 5
    smoothing = 0;
end
if nargin < 4
    % default 1d-statistic:
    stat1d = 1;
end
if nargin < 3
    % default amount of projection lines:
    directions = 100;
    if dim == 1
        directions = 1;
    end
end

% just in case, check the amount of directions
if directions > max_directions
    disp(['Too many directions! Reduced the directions to ',...
        num2str(max_directions),'.']);
    directions = max_directions;
end

% compute the statistic!
[dn_sum, dn_max] = dt_gof_statistic(samples_exp, samples_model, [], [], ...
    directions, stat1d, two_sample, smoothing, use_normalised_stat);

% the average of the discrepancies to different directions:
D_sum = d1_scaling * dn_sum; 
% the maximum discrepancy of the different directions:
D_max = d2_scaling * dn_max; 

end




