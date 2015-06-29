function[] = brainmap();

thr = 0%1.6; % Z threshold
thr_num_exp_bd = 0;
thr_num_exp_pc = 0;
sort_by_mean = 1;
fd = load('FunctionalDecoding.mat');

h = figure('units','normalized','outerposition',[0 0 1 1]);

% BD
subplot(1,2,1);
bd = (fd.BDprofile_LR .* (fd.BDprofile_Z > thr))';
bd = bd(:,1:17);
bd = bd(fd.BD_Nexperiments > thr_num_exp_bd,:);
notzero = find(sum(bd,2));
bd = bd(notzero,:);

if sort_by_mean == 1
    % For sorting by mean:
    [i,j] = find(bd);
    clear m
    for k = 1:length(bd)
        m(k) = mean(j(i == k));
    end
    [s,ind] = sort(m);
else 
    % For sorting by Fiedler:
    m = laplacian_eigen(bd);  
    [s,ind] = sort(m(:,1));
end
ind = ind(isfinite(s));
[~,maxi] = max(mean(bd(ind,1),2));
if maxi > 5
    ind = fliplr(ind);
end

imagesc(bd(ind,:));
hold on;
set(gca,'YTick',[1:length(ind)],'YTickLabel', fd.BDnames(notzero(ind)));
set(gca,'XTick',[1:17],'XTickLabel',num2cell([5:5:85]));
title('behavioral domains (BD)');
colorbar;

% PC
subplot(1,2,2);
pc = (fd.PCprofile_LR .* (fd.PCprofile_Z > thr))';
pc = pc(:,1:17);
pc = pc(fd.PC_Nexperiments > thr_num_exp_pc,:);
notzero = find(sum(pc,2));
pc = pc(notzero,:);

if sort_by_mean == 1
    % For sorting by mean:
    [i,j] = find(pc);
    clear m
    for k = 1:length(pc)
        m(k) = mean(j(i == k));
    end
    [s,ind] = sort(m);
else
    % For sorting by Fiedler:
    m = laplacian_eigen(pc);    
    [s,ind] = sort(m(:,1));
end
ind = ind(isfinite(s));
[~,maxi] = max(mean(pc(ind,1),2));
if maxi > 5
    ind = fliplr(ind);
end

imagesc(pc(ind,:));
hold on;
set(gca,'YTick',[1:length(ind)],'YTickLabel', fd.PCnames(notzero(ind)));
set(gca,'XTick',[1:17],'XTickLabel',num2cell([5:5:85]));
title('PC');
colorbar;

saveas(h,'fig.brainmap.pdf','pdf')

