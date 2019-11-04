load('nback_remOutliers.mat'); %load data

aut_a = groups(aut,aut_a); %run this code on autism participants
con_a = groups(con,con_a); %run this code on controls

save('nback_rats.mat','aut_a','con_a')

function [sub_a] = groups(sub,sub_a)

    subRatHit = nan(length(sub),1); %ratios of distractor to no distractor amplitudes (hits) (initialisation of variable)
    subRatFA = nan(length(sub),1); %ratios of distractor to no distractor amplitudes (FAs) (initialisation of variable)
    subRatMiss = nan(length(sub),1); %ratios of distractor to no distractor maplitudes (misses) (initialisation of variable)
    for n = 1:length(sub)

       %for each participant
       sub_a(n).RatHit = sub_a(n).taskIrf_hit_mod/sub(n).taskIrf_hit_mod;
       sub_a(n).RatFA = sub_a(n).taskIrf_FA_mod/sub(n).taskIrf_FA_mod;
       sub_a(n).RatMiss = sub_a(n).taskIrf_miss_mod/sub(n).taskIrf_miss_mod;

    end
    
end