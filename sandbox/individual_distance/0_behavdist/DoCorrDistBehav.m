function [] = DoCorrDistBehav()

load data_dist
addpath(genpath('../../../utils'));
addpath(genpath('../../../../00_matlab_packages/palm-alpha52/'));

%% variables:
raw = dlmread('subs_2.csv');
sublist = raw(:,2);
cols = [47:70 75:78];
covs = raw(:,cols);
cov_names = {   'Mean_thresh_audition',         'Mean_thresh_vision',...
                'Mean_taskfocus_audition',      'Mean_taskfocus_vision',...
                'Mean_future_audition',         'Mean_future_vision', ...
                'Mean_past_audition',           'Mean_past_vision', ...
                'Mean_self_audition',           'Mean_self_vision',...
                'Mean_other_audition',          'Mean_other_vision', ...
                'Mean_NegPos_audition',         'Mean_NegPos_vision', ...
                'Mean_images_audition',         'Mean_images_vision',...
                'Mean_words_audition',          'Mean_words_vision',...
                'Mean_intrusive_audition',      'Mean_intrusive_vision',...
                'Mean_SpecVag_audition',        'Mean_SpecVag_vision',...
                'Mean_thresh_audition_zscore',  'Mean_thresh_vision_zscore',...
                'Mean_dprime_audition',         'Mean_dprime_vision',...
                'age',                          'sex'};

%% prepare distance data
ind = 1:size(data,1);
d = data(ind,:,:,:);


% check for missing datapoints
ind = (find(min([sum(squeeze(data(:,1,1,:))')' sum(squeeze(data(:,2,1,:))')' sum(squeeze(data(:,1,2,:))')' sum(squeeze(data(:,2,2,:))')']')'));
d = data(ind,:,:,:);
sublist = sublist(ind);
covs = raw(ind,cols);

% normalize so peak value is 1
for i = 1:length(ind)
    for h = 1:2
        for l = 1:2
            dn(i,h,l,:) = d(i,h,l,:) ./ max(d(i,h,l,:));
        end
    end
end

surf = SurfStatReadSurf({['../../../../yeoTopo/lh.pial'],['../../../../yeoTopo/rh.pial']});



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

    Mean_thresh_audition = term(covs(:,1));
    Mean_thresh_vision = term(covs(:,2));
    Mean_taskfocus_audition = term(covs(:,3));
    Mean_taskfocus_vision = term(covs(:,4));
    Mean_future_audition = term(covs(:,5));
    Mean_future_vision = term(covs(:,6));
    Mean_past_audition = term(covs(:,7));
    Mean_past_vision = term(covs(:,8));
    Mean_self_audition = term(covs(:,9));
    Mean_self_vision = term(covs(:,10));
    Mean_other_audition = term(covs(:,11));
    Mean_other_vision = term(covs(:,12));
    Mean_NegPos_audition = term(covs(:,13));
    Mean_NegPos_vision = term(covs(:,14));
    Mean_images_audition = term(covs(:,15));
    Mean_images_vision = term(covs(:,16));
    Mean_words_audition = term(covs(:,17));
    Mean_words_vision = term(covs(:,18));
    Mean_intrusive_audition = term(covs(:,19));
    Mean_intrusive_vision = term(covs(:,20));
    Mean_SpecVag_audition = term(covs(:,21));
    Mean_SpecVag_vision = term(covs(:,22));
    Mean_thresh_audition_zscore = term(covs(:,23));
    Mean_thresh_vision_zscore = term(covs(:,24));
    Mean_dprime_audition = term(covs(:,25));
    Mean_dprime_vision = term(covs(:,26));
    age = term(covs(:,27));
    sex = term(covs(:,28));

ds = SurfStatSmooth(d, surf, 10); % smoothing of 20 mm
total_sa = sum(ds,2);

% To get contrasts of seed points:
for i = 1:size(covs,2)    
    c = covs(:,i);    
    % M = 1 + term(total_sa) + term(c);    
    M = 1 + term(c);    
    for j = 1:2
        input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
        %input = ds;            
        slm  = SurfStatLinMod(input, M, surf);
        contrast = M(2);
        %contrast = M(3);
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


for i = [1:size(covs,2)]
    input = [logical(stats1{i}.clusid) .* logical(stats2{i}.clusid) - logical(stats_n1{i}.clusid) .* logical(stats_n2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_dist_same_contrast_' cov_names{i} '.png']);
        close all;   
    end
end
for i = [1:size(covs,2)]
    input = [logical(stats1{i}.clusid) .* logical(stats_n2{i}.clusid) - logical(stats_n1{i}.clusid) .* logical(stats2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_dist_Other_same_contrast_' cov_names{i} '.png']);
        close all;   
    end
end
for i = 1:2:size(covs,2)-1
    disp(i);
    input = [logical(stats1{i+1}.clusid) .* logical(stats2{i}.clusid) - logical(stats_n1{i+1}.clusid) .* logical(stats_n2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_dist_modality_specific_' cov_names{i} '.png']);
        close all;
    end
end
for i = 1:2:size(covs,2)-1
    input = [logical(stats1{i+1}.clusid) .* logical(stats_n2{i}.clusid) - logical(stats_n1{i+1}.clusid) .* logical(stats2{i}.clusid)];
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_dist_Other_modality_specific_' cov_names{i} '.png']);
        close all;
    end
end


map = {'vis','aud'};
for i = 1:2:size(covs,2)-1 
        input = [logical(stats2{i}.clusid) - logical(stats_n2{i}.clusid)];
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_contrast_' cov_names{i} '_map_' map{2} '.png']);
        close all;   
end
for i = 2:2:size(covs,2)
        input = [logical(stats1{i}.clusid) - logical(stats_n1{i}.clusid)];
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_contrast_' cov_names{i} '_map_' map{1} '.png']);
        close all;   
end




for i = 3:4 %[1:size(covs,2)]
    input = [logical(stats1{i}.clusid) - logical(stats_n1{i}.clusid) ];
    %input = [logical(stats_n1{i}.clusid) ] .* -1;
    if sum(input) ~= 0 
        h = figure('visible','off'); 
        SurfStatView(input, surf, cov_names{i});
        SurfStatColormap(colormap(rgb));
        SurfStatColLim([-1 1]);
        saveas(gcf,['Fig_SA_same_contrast_' cov_names{i} '.png']);
        close all;   
    end
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

%% 
i = 4; j = 1; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm4v = slm;

i = 3; j = 1; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm3v = slm;

i = 4; j = 2; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm4a = slm;

i = 3; j = 2; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm3a = slm;

figure; SurfStatView(slm4v.t-slm4a.t, surf);
saveas(gcf,'Fig.unthresh.slm4v-slm4a.png');
SurfStatView(slm3v.t-slm3a.t, surf);
saveas(gcf,'Fig.unthresh.slm3v-slm3a.png');
SurfStatView(mean([slm3v.t slm4v.t],1)-mean([slm3a.t slm4a.t],1), surf);
saveas(gcf,'Fig.unthresh.slm_v-slm_a.png');
SurfStatView(slm4v.t-slm3a.t, surf);
saveas(gcf,'Fig.unthresh.slm_v-slm_a.modspecific.png');

figure; SurfStatView(slm4v.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm4v.png');
SurfStatView(slm3a.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm3a.png');
SurfStatView(slm4a.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm4a.png');
SurfStatView(slm3v.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm3v.png');



i = 4; j = 1; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm4v = slm;

i = 3; j = 1; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm3v = slm;

i = 4; j = 2; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm4a = slm;

i = 3; j = 2; c = covs(:,i); M = 1 + term(c);
input = squeeze([squeeze(d(:,1,j,:)) squeeze(d(:,2,j,:))]);            
slm  = SurfStatLinMod(input, M, surf); contrast = M(2);
slm  = SurfStatT(slm, contrast); 
slm3a = slm;

figure; SurfStatView(slm4v.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm4v.png');
SurfStatView(slm3a.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm3a.png');
SurfStatView(slm4a.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm4a.png');
SurfStatView(slm3v.t, surf)
saveas(gcf,'Fig.hcp.unthresh.slm3v.png');