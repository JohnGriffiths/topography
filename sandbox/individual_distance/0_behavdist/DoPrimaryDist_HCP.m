function [] = DoPrimaryDist_HCP(sub)

run = 'dist';
dir = '/a/documents/connectome/_all/'; 

labels = [46 34]; % v1: S_calcarine ; a1: G_temp_sup-G_T_transv
hemi = {'L','R'};

addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/surfstat'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/geodesic'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils'));
addpath(genpath('/scr/animals1/Dropbox/01_code/00_matlab_packages/palm-alpha52/'));

switch run
    case 'dist'
        filename = ['/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/dist.hcp.' sub '.v1.a1.mat'];
        if ~exist(filename, 'file');
            disp(sub);
            distances = [];
            for h = 1:2
                fsa_reg = readGiiSurf([dir sub '/MNINonLinear/fsaverage_LR32k/' sub '.' hemi{h} '.sphere.32k_fs_LR.surf.gii']);                       
                sub_reg = readGiiSurf([dir sub '/MNINonLinear/Native/' sub '.' hemi{h} '.sphere.reg.reg_LR.native.surf.gii']);      
                indVal = zeros(1,length(fsa_reg.coord));
                for i = 1:length(fsa_reg.coord)   
                    [~, indVal(i)] = min(sum(bsxfun(@minus,fsa_reg.coord(:,i), sub_reg.coord).^2).^0.5);                     
                end
                fslabel = gifti([dir sub '/MNINonLinear/Native/' sub '.' hemi{h} '.aparc.a2009s.native.label.gii']);
                for j = 1:length(labels)
                    prim = find(fslabel.cdata == fslabel.labels.key(labels(j)));        
                    dist = distExactGeodesic(prim, 'Native', hemi{h}, 'distance', dir , sub);
                    distances(h,j,:) = dist(indVal);
                end                             
            end    
            disp(filename);
            save(filename,'-v7.3','distances');
        end


    case 'data'
        raw = csvread('../../../../../02_data/hcp/covs/sustained_attention_task.csv',1,0);     
        sublist = raw(:,1);
        
        data = zeros(length(sublist),2,2,32492);
        count = 1;
        for s = 1:length(sublist)
            filename = ['dist.hcp.' num2str(sublist(s)) '.v1.a1.mat'];
            if exist(filename,'file');
                d = load(filename);
                data(s, :, :, :) = d.distances;
                disp(count); count = count + 1;
            end
        end
        save('data_dist_hcp.mat','-v7.3','data');
end
