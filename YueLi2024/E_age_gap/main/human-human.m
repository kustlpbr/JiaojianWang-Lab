clc;clear;
%input data - human data
all_mats_human = load('human_wm_gm.mat');
all_mats_human = struct2array(all_mats)';
all_age_human = load('8_12_age.mat');
all_age_human = struct2array(all_age_human);
thresh = 0.01; % 0.05
no_sub_human = 370; % the number of subjects
no_node1_human = 260; % the number of parameters
no_node2_human = 1;
%train and building model
%10-cross-validation
for j=1:100;
    fprintf('\n subj #%6.3f',j);
    x = 1:370;
    n = 333;
    random_num(j,:) = x(randperm(numel(x),n));
    random_num(j,:) = sort(random_num(j,:));
    %random_num = random_num';
    for i=1:333;
        a = random_num(j,i);
        train_mats(:,:,i) = all_mats(:,:,a);
        train_age(i,:) = all_age(a,:);
    end
    train_mats = reshape(train_mats,[260,1,333]);
    train_vcts = reshape(train_mats,[260,333]);
    no_sub = size(train_mats,3);
    no_node1 = size(train_mats,1);
    no_node2 = size(train_mats,2);
    %correlate
    [r_mat,p_mat] = corr(train_vcts',train_age);
    r_mat = reshape(r_mat,no_node1,no_node2);
    p_mat = reshape(p_mat,no_node1,no_node2);
    pos_mask =zeros(no_node1,no_node2);
    neg_mask =zeros(no_node1,no_node2);
    pos_edges = find(r_mat > 0.05 & p_mat < thresh);
    neg_edges = find(r_mat < -0.05 & p_mat < thresh);
    pos_mask(pos_edges) = 1;
    neg_mask(neg_edges) = 1;
    %get sum of all edge in train subs
    %(divide by 2 to control for the fact that matrics are symmetric)
    train_sumpos = zeros(no_sub,1);
    train_sumneg = zeros(no_sub,1);
    for ss = 1:size(train_sumpos);
        train_sumpos(ss) = sum( sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum( sum(train_mats(:,:,ss).*neg_mask))/2;
    end
    %build model on train subs
    b = regress(train_age , [train_sumpos,train_sumneg , ones(no_sub,1)]);
    %run model on test sub
    for i=1:no_sub_human;
        fprintf('\n leaving out subj #%6.3f',i);
        test_mat = all_mats_human(:,:,i);
        test_sumpos = sum(sum(test_mat.*pos_mask))/2;
        test_sumneg = sum(sum(test_mat.*neg_mask))/2;
        behav_pred(i,1) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    end
    all_behav_pred(:,j) = behav_pred(:,1);
end
all_behav_pred = all_behav_pred';
pre = average(all_behav_pred)';
[R,P] = corr(pre,all_age_human);
MAE = mae(pre,all_age_human);
%predict data - actual data
subtraction = abs(pre - all_age_human);
[subtraction_R,subtraction_P] = corr(subtraction,all_age_human)
MAE_subtraction = mae(subtraction,all_age_human)

figure(1);plot(all_age_human,pre,'r.');lsline;
figure(2);plot(all_age_human,subtraction,'r.');lsline;

