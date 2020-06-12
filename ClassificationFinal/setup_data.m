
clearvars;
load('../data/nback_remOutliers.mat')

n_tp = length(aut(1).taskIrf_miss);

for distractors = {'', '_a'}
    for dtype = {'hit', 'FA', 'miss'}

        if length(distractors{1}) == 0
            fprintf('\nno distractors\n')
            aut_dat = aut;
            con_dat = con;
        else
            fprintf('\ndistractors\n')
            aut_dat = aut_a;
            con_dat = con_a;
        end

        X = [];
        Y = [];
        subnums = [];
        
        for sub = 1:length(aut_dat)
            if ((aut(sub).hitRate == 1) || (aut_a(sub).hitRate == 1)) && strcmp(dtype{1}, 'miss')
                fprintf('\n bad sub-%02d for miss type 1', sub)
                continue
            end
            if ((aut(sub).faRate == 0) || (aut_a(sub).faRate == 0)) && strcmp(dtype{1}, 'FA')
                fprintf('\n bad sub sub-%02d for FA  type 1')
                continue
            end
            if any(isnan(aut(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(aut_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
                fprintf('\n bad sub-%02d for %s type 1', dtype{1})
                continue
            end
            X = [X; aut_dat(sub).(sprintf('taskIrf_%s', dtype{1}))];
            Y = [Y; 1];
            subnums=[subnums; sub];
        end
        for sub = 1:length(con_dat)
            if ((con(sub).hitRate == 1) || (con_a(sub).hitRate == 1)) && strcmp(dtype{1}, 'miss')
                fprintf('\n bad sub for miss type 0')
                continue
            end
            if ((con(sub).faRate == 0) || (con_a(sub).faRate == 0)) && strcmp(dtype{1}, 'FA')
                fprintf('\n bad sub for FA  type 0')
                continue
            end
            if any(isnan(con(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(con_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
                fprintf('\n bad sub for %s type 0', dtype{1})
                continue
            end
            X = [X; con_dat(sub).(sprintf('taskIrf_%s', dtype{1}))];
            Y = [Y; 0];
            subnums=[subnums; sub];
        end
        data = struct('X',X,'Y',Y, 'subnum',subnums);
        save(sprintf('XY%s_%s.mat', distractors{1}, dtype{1}), 'data')
    end
end

%% take difference
for dtype = {'hit', 'FA', 'miss'}
    
    % no distractor (not a)
    X = [];
    Y = [];
    subnums=[];

    for sub = 1:length(aut)
        if ((aut(sub).hitRate == 1) || (aut_a(sub).hitRate == 1)) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub-%02d for miss type 1', sub)
            continue
        end
        if ((aut(sub).faRate == 0) || (aut_a(sub).faRate == 0)) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub-%02d for FA  type 1', sub)
            continue
        end
        if any(isnan(aut(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(aut_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub-%02d for %s type 1 (NaNs)', sub, dtype{1})
            continue
        end
        dat = mad(aut(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        dat_a = mad(aut_a(sub).(sprintf('taskIrf_%s', dtype{1})),1);
%         X = [X; (dat-dat_a)./(dat+dat_a)];
        X = [X; dat_a-dat];
        Y = [Y; 1];
    end
    for sub = 1:length(con)
        if ((con(sub).hitRate == 1) || (con_a(sub).hitRate == 1)) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub-%02d for miss type 0', sub)
            continue
        end
        if ((con(sub).faRate == 0) || (con_a(sub).faRate == 0)) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub-%02d for FA  type 0', sub)
            continue
        end
        if any(isnan(con(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(con_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub-%02d for %s type 0 (NaNs)', sub, dtype{1})
            continue
        end
        dat = mad(con(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        dat_a = mad(con_a(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        X = [X; dat_a-dat];
%         X = [X; (dat-dat_a)./(dat+dat_a)];
        Y = [Y; 0];
        
    end
    data = struct('X',X,'Y',Y, 'subnum',subnums);
    save(sprintf('XY_distractors-diff_MAD_%s.mat', dtype{1}), 'data')
end


%% old code for trying temporal decoding in matlab (moved to Python, see classify.py)

%{
bin_width = 10;
accs = zeros(20,70);
for iter = 1:20
    for t0 = 1:70
        SVMModel = fitcsvm(X(:,t0:t0+bin_width),Y,'Standardize',true);
        CVSVMModel = crossval(SVMModel);
        accs(iter, t0) = 1-kfoldLoss(CVSVMModel, 'loss', 'classiferr');
    end
end

figure; 
subplot(2,1,1)
hold on
errorbar(mean(X(1:22,:),1), std(X(1:22,:),[],1)/sqrt(22)); errorbar(mean(X(23:end,:),1), std(X(23:end,:),[],1)/sqrt(24))
xlabel('Time point')
ylabel('Mean IRF')
subplot(2,1,2)
hold on
errorbar(mean(accs,1), std(accs,[],1)/sqrt(20))
plot([0,80], [.5, .5], 'r')
xlabel(sprintf('Start of %d time point window', bin_width))
ylabel('SVM accuracy')
%}
