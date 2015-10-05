function [] = DoPrimaryDist(sub)

run = 'dist';

%dir='/scr/tantalum1/freesurfer_backup/';
dir = '/scr/animals1/lsd_freesurfer/'; %/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/';

%toRun = [];
%    for i = 1:length(toRun)
%      filename = [dir num2str(toRun(i)) '/surf/lh.pial'];
%      if exist(filename,'file');
%         DoPrimaryDist(num2str(toRun(i)));
%      end
%    end
%raw = dlmread('subs_2.csv');
%sublist = raw(:,2);

labels = [46 34]; % v1: S_calcarine ; a1: G_temp_sup-G_T_transv
hemi = {'lh','rh'};

addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/surfstat'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/geodesic'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils'));

switch run
    case 'dist'
        % for s = 1:length(sublist)
            % sub = num2str(sublist(s));
            % sub = sublist{s};
            filename = ['/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/dist.' sub '.v1.a1.mat'];
            if ~exist(filename, 'file');% && exist([dir sub '/surf/lh.sphere.reg'], 'file') ;
                disp(sub);
                distances = [];
                for h = 1:2
                    fsa_reg = SurfStatReadSurf(['/scr/tantalum1/freesurfer_backup/fsaverage5/surf/' hemi{h} '.sphere.reg']);       
                    disp([dir sub '/surf/' hemi{h} '.sphere.reg']);
                    sub_reg = SurfStatReadSurf([dir sub '/surf/' hemi{h} '.sphere.reg']);        
                    indVal = zeros(1,length(fsa_reg.coord));
                    for i = 1:length(fsa_reg.coord)   
                        [~, indVal(i)] = min(sum(bsxfun(@minus,fsa_reg.coord(:,i), sub_reg.coord).^2).^0.5);            
                    end
                    [~, fslabel, colortable] = read_annotation([ dir sub '/label/' hemi{h} '.aparc.a2009s.annot']);
                    for j = 1:length(labels)
                        prim = find(fslabel == colortable.table(labels(j),5));        
                        dist = distExactGeodesic(prim, 'freesurfer', hemi{h}, 'distance', dir , sub);
                        distances(h,j,:) = dist(indVal);
                    end                             
                end    
                disp(filename);
                save(filename,'-v7.3','distances');
            end
        % end

    case 'data'
        data = zeros(length(sublist),2,2,10242);
        count = 1;
        for s = 1:length(sublist)
            filename = ['dist.' num2str(sublist(s)) '.v1.a1.mat'];
            if exist(filename,'file');
                d = load(filename);
                data(s, :, :, :) = d.distances;
                %disp(count); 
                count = count + 1;
            else 
                disp(sublist(s));
            end
        end
        disp(num2str(count - 1));
        save('data_dist.mat','-v7.3','data');
end

