function [data] = DoGetDistances()

dataDir = '/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands';

raw = dlmread('subs.csv','\t');
sublist = raw(:,2);

hemi = {'lh', 'rh'};
label = {'S_calcarine_fsa5','G_temp_sup-G_T_transv_fsa5'};
data = zeros(length(sublist),2,2,10242);

for s = 1:length(sublist)
    sub = num2str(sublist(s));
    for h = 1:length(hemi)
        filename = [dataDir '/' sub '/distance_maps/' sub '_' hemi{h} '_geoDist_fsa5.mat'];
        if exist(filename,'file');
            d = load(filename);
            dist = d.dataAll;

            for l = 1:length(label)
                fname = [dataDir '/' sub '/labels/fsa5/' hemi{h} '.' label{l} '.label'];

                fid = fopen(fname, 'r');
                fgets(fid) ;       
                line = fgets(fid) ;
                nv = sscanf(line, '%d') ;
                lraw = fscanf(fid, '%d %f %f %f %f\n') ;
                lraw = reshape(lraw, 5, nv) ;
                lraw = lraw' + 1;
                fclose(fid) ;

                data(s, h, l, :) = min(dist(lraw(:,1),:));
            end
        end
    end
    disp(sub);
end

save('data_dist.mat','-v7.3','data');