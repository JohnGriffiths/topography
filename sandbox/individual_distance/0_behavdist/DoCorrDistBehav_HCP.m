function [] = DoCorrDistBehav_HCP()

load data_dist_hcp
addpath(genpath('../../../utils'));
addpath(genpath('../../../../00_matlab_packages/palm-alpha52/'));

%% variables:
raw = csvread('../../../../../02_data/hcp/covs/sustained_attention_task.csv',1,0);     
sublist = raw(:,1);

covs = raw(:,2:end);
cov_names = {'FS_InterCranial_Vol' 'SCPT_TP' 'SCPT_TN' 'SCPT_FP' 'SCPT_FN' 'SCPT_TPRT' 'SCPT_SEN' 'SCPT_SPEC' 'SCPT_LRNR' 'Flanker_Unadj' 'Flanker_AgeAdj' 'CarSort_Unadj' 'CardSort_AgeAdj'};

%% prepare distance data

% check for missing datapoints
ind = (find(min([sum(squeeze(data(:,1,1,:))')' sum(squeeze(data(:,2,1,:))')' sum(squeeze(data(:,1,2,:))')' sum(squeeze(data(:,2,2,:))')']')'));
ind = ind(find(covs(ind,2) ~= 0));
d = data(ind,:,:,:);
sublist = sublist(ind);
covs = raw(ind,2:end);

% normalize so peak value is 1
% for i = 1:length(ind)
%     for h = 1:2
%         for l = 1:2
%             dn(i,h,l,:) = d(i,h,l,:) ./ max(d(i,h,l,:));
%         end
%     end
% end

surfL = gifti('../../../data/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii');
surfR = gifti('../../../data/Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii');
surf.coord = [surfL.vertices; surfR.vertices]';
surf.tri = [surfL.faces; surfR.faces + 32492];


%% run regressions:
% can also run with 'dn', which is normalized

% sub = 100307;
% aparc = gifti(['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
% '.L.aparc.a2009s.32k_fs_LR.label.gii']);
% aparc = aparc.cdata;
% noncortex = zeros([1 length(surf1.vertices)]);
% noncortex(find(aparc == 0)) = 1;
% cortex = zeros([1 length(surf1.vertices)]);
% cortex(find(noncortex == 0)) = 1;
% mask = cortex;

    FS_InterCranial_Vol = term(covs(:,1));
    SCPT_TP = term(covs(:,2));
    SCPT_TN = term(covs(:,3));
    SCPT_FP = term(covs(:,4));
    SCPT_FN = term(covs(:,5));
    SCPT_TPRT = term(covs(:,6));
    SCPT_SEN = term(covs(:,7));
    SCPT_SPEC = term(covs(:,8));
    SCPT_LRNR = term(covs(:,9));   
    Flanker_Unadj = term(covs(:,10));
    Flanker_AgeAdj = term(covs(:,11));
    CardSort_Unadj = term(covs(:,12));
    CardSort_AgeAdj = term(covs(:,13)); 
    
% To get contrasts of seed points:
for i = 1:size(covs,2)    
    c = covs(:,i+1)-1;
    M = 1 + FS_InterCranial_Vol + term(c);    
    for j = 1:2
        input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
        slm  = SurfStatLinMod(input, M, surf);
        contrast = M(3);
        slm  = SurfStatT(slm, contrast); % Take this + and -  to look at both directions!
        %mask = ones(size(input,1),1);
        %slm.t= slm.t.*mask; %tvalue    
        mask = ones(size(input,2),1);
        clusp = 0.05 ;        
        if j == 1
            [stats1{i}.pval,stats1{i}.peak,stats1{i}.clus,stats1{i}.clusid] = SurfStatP(slm, mask, clusp);
        else
            [stats2{i}.pval,stats2{i}.peak,stats2{i}.clus,stats2{i}.clusid] = SurfStatP(slm, mask, clusp);
        end      
        contrast = contrast * -1;
        slm  = SurfStatT(slm, contrast);
        if j == 1
            [stats_n1{i}.pval,stats_n1{i}.peak,stats_n1{i}.clus,stats_n1{i}.clusid] = SurfStatP(slm, mask, clusp);
        else
            [stats_n2{i}.pval,stats_n2{i}.peak,stats_n2{i}.clus,stats_n2{i}.clusid] = SurfStatP(slm, mask, clusp);
        end  
        clear slm
    end
    disp(i);
end

%% make figures: 
% check contrast 5:
rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

for i = [1:size(covs,2)-1]
    input = [logical(stats1{i}.clusid) .* logical(stats2{i}.clusid) - logical(stats_n1{i}.clusid) .* logical(stats_n2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i+1});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig.hcp.inlcd_Vol.2.' cov_names{i} '.png']);
        close all;   
    end
end
for i = [1:size(covs,2)-1]
    input = [logical(stats1{i}.clusid) .* logical(stats_n2{i}.clusid) - logical(stats_n1{i}.clusid) .* logical(stats2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig.hcp.other.inlcd_Vol.2.' cov_names{i+1} '.png']);
        close all;   
    end
end

map = {'vis','aud'};
for i = 1:size(covs,2)-1
        input = [logical(stats2{i}.clusid) - logical(stats_n2{i}.clusid)];
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i+1});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig.hcp.inlcd_Vol.2.' cov_names{i+1} '_map_' map{2} '.png']);
        close all;   
end
for i = 1:size(covs,2)-1
        input = [logical(stats1{i}.clusid) - logical(stats_n1{i}.clusid)];
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i+1});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig.hcp.inlcd_Vol.2.' cov_names{i+1} '_map_' map{1} '.png']);
        close all;   
end



% 
% %%
% %% write out file:
% cov_header(1) = cellstr(['subID']);
% for i = 1:length(cov_names)
%     cov_header(1 + i) = cellstr(cov_names(i));
% end
% for i = 1:20484
%     cov_header(25 + i) = cellstr(['vis_' num2str(i)]);
% end
% for i = 1:20484
%     cov_header(25 + 20484 + i) = cellstr(['aud_' num2str(i)]);
% end
% 
% data_all(:,1) = sublist;
% data_all(:,2:25) = covs;
% data_all(:,26:26+length(dbv)-1) = dbv;
% data_all(:,26+length(dbv):26+length(dbv)+20484-1) = dbv;
% 
% csvwriteh('dist_data_all.csv', data_all, cov_header);
% 
% %% Regression
% for j = 1:2:24
%     for i = 1:20484
%         [bv(j/2+0.5,i,:),~,~,~,sv(j/2+0.5,i,:)] = regress(dbv(:,i),[ones(length(covs),1) covs(:,j) covs(:,j+1)]);
%         [ba(j/2+0.5,i,:),~,~,~,sa(j/2+0.5,i,:)] = regress(dba(:,i),[ones(length(covs),1) covs(:,j) covs(:,j+1)]);        
%     end
%     disp(j);
% end

% dbv = [squeeze(d(:,1,1,:)) squeeze(d(:,2,1,:))];
% dba = [squeeze(d(:,1,2,:)) squeeze(d(:,2,2,:))];
% rv = corr(dbv, covs);
% ra = corr(dba, covs);
% rv(isnan(rv)) = 0;
% ra(isnan(ra)) = 0;
% 
% a = pca(rv);
% figure; SurfStatView(a(:,1),surf);
% figure; SurfStatView(a(:,2),surf);
% b = pca(ra);
% figure; SurfStatView(b(:,1),surf);
% figure; SurfStatView(b(:,2),surf);
% 
% 
% dnbv = [squeeze(dn(:,1,1,:)) squeeze(dn(:,2,1,:))];
% dnba = [squeeze(dn(:,1,2,:)) squeeze(dn(:,2,2,:))];
% rnv = corr(dnbv, covs);
% rna = corr(dnba, covs);
% rnv(isnan(rnv)) = 0;
% rna(isnan(rna)) = 0;
% 
% a = pca(rnv);
% figure; SurfStatView(a(:,1),surf);
% figure; SurfStatView(a(:,2),surf);
% b = pca(rna);
% figure; SurfStatView(b(:,1),surf);
% figure; SurfStatView(b(:,2),surf);
% 
% figure;
% for i = 1:24
%     SurfStatView(ra(:,i),surf, cov_names{i});
%     pause(5);
% end

i = 3;
input = [logical(stats1{i+1}.clusid) .* logical(stats2{i}.clusid) - logical(stats_n1{i+1}.clusid) .* logical(stats_n2{i}.clusid)];
input2 = squeeze([squeeze(d(:,1,1,:)) squeeze(d(:,2,1,:))]);
figure; scatter(mean(input2(:,find(input)),2), covs(:,4), 'rx');
input2 = squeeze([squeeze(d(:,1,2,:)) squeeze(d(:,2,2,:))]);
hold on; scatter(mean(input2(:,find(input)),2), covs(:,3), 'bx');


