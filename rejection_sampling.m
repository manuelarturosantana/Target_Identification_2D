
load("Identification/classification_results/plane_results_averaged_tks.mat")
inds_miss = results.trials.correct == 0;

%% TODO: This is difficult. It seems that many of the good classifications are near the bad ones.
%%
%%
pred_idx = results.trials.pred_idx;
lin_idx = sub2ind(size(results.reject), (1:numel(pred_idx))', pred_idx);
rej = results.reject(lin_idx);

n_reject_correct = sum(rej(inds_miss==1));
n_reject_false   = sum(rej(inds_miss==0));

n_reject_correct / sum(inds_miss)
n_reject_false / sum(~inds_miss)

%%
n_reject_correct  = 0;
maxes = [];
smaxes = [];
for ii = 1:length(inds_miss)
        s = results.scores(ii,:);
        m1 = max(s); maxes(end+1) = m1;
        m2 = max(s(s~=m1)); smaxes(end+1) = m2;
        thres = mean(m1 ./ s(s~=m1)) + var(m1 ./ s(s~=m1));
   
        if (m1 / m2) < thres
            n_reject_correct = n_reject_correct + 1;
        end
end


n_reject_correct / sum(inds_miss)
%% What if we just try if m2 is within a certain percentage of m1
percent = 0.2;
n_reject_correct  = 0;
n_reject_false   = 0;

for ii = 1:length(inds_miss)
        s = results.scores(ii,:);
        m1 = max(s); 
        m2 = max(s(s~=m1)); 
        
   
        if ((m1 - m2) / m1) > 0.02
            if (inds_miss(ii) == 0)
                n_reject_correct = n_reject_correct + 1;
            else 
                n_reject_false = n_reject_false + 1;
            end
        end

end


n_reject_correct / sum(inds_miss)
n_reject_false / sum(~inds_miss)

%%
figure(1)
clf
plot(mean(results.scores,2).')
figure(2)
clf
plot(var(results.scores, 0, 2).')

%%
rel_diff_correct = [];
rel_diff_wrong = [];
for ii = 1:length(inds_miss)
    s = results.scores(ii,:);
    m1 = max(s); 
    m2 = max(s(s~=m1));

    if inds_miss(ii) == 0
    
        rel_diff_correct(end+1) = abs(m1 - m2) / m1;
    else
         rel_diff_wrong(end+1) = abs(m1 - m2) / m1;
    end
end

mean(rel_diff_wrong)
mean(rel_diff_correct)
var(rel_diff_wrong)
var(rel_diff_correct)

%% Overall diagonises sofar. This is tough. it doesn't appear like there is a good
% easy rule to get rid of false positives. I think chat gpt is right
% looking at the statistics more...