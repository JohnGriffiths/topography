function[] = DoSurfArea(sub)

addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/surfstat'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils/geodesic'));
addpath(genpath('/scr/animals1/Dropbox/01_code/topography/utils'));
addpath(genpath('/scr/animals1/Dropbox/01_code/00_matlab_packages/palm-alpha52/'));

% prepare surfarea data
dir = '/scr/animals1/lsd_freesurfer/'; 
hemi = {'lh','rh'};

filename = ['/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/surfarea.' sub '.mat'];
if ~exist(filename, 'file');
    sa = [];
    for h = 1:2        
        rawg = read_curv([dir sub '/surf/' hemi{h} '.area.pial']);

        fsa_reg = SurfStatReadSurf(['/scr/tantalum1/freesurfer_backup/fsaverage5/surf/' hemi{h} '.sphere.reg']);       
        disp([dir sub '/surf/' hemi{h} '.sphere.reg']);
        sub_reg = SurfStatReadSurf([dir sub '/surf/' hemi{h} '.sphere.reg']);        
        indVal = zeros(1,length(fsa_reg.coord));
        for i = 1:length(fsa_reg.coord)   
            [~, indVal(i)] = min(sum(bsxfun(@minus,fsa_reg.coord(:,i), sub_reg.coord).^2).^0.5);        
        end

        [~, zones] = distExactGeodesic(indVal, 'freesurfer_pial', hemi{h}, 'zones', dir, sub);
        zonesAll = zeros(length(zones),1);
        zonesAll(indVal) = [1:length(indVal)];
        for i = 1:length(indVal)
            zonesAll(find(zones == zones(indVal(i)))) = i;
        end    
        satot = zeros(length(fsa_reg.coord),1); 
        for i = 1:length(fsa_reg.coord) 
            satot(i) = sum(rawg(find(zonesAll == i)));        
        end
        satot(find(satot > 150)) = 0;
        sa = [sa; satot];
    end
    save(filename,'-v7.3','sa');
end

%%
% data = zeros(length(sublist),10242*2);
% count = 1;
%         for s = 1:length(sublist)
%             filename = ['surfarea.' num2str(sublist(s)) '.mat'];
%             if exist(filename,'file');
%                 d = load(filename);
%                 data(s, :, :, :) = d.sa;
%                 disp(count); count = count + 1;
%             end
%         end
%         save('data_surfarea.mat','-v7.3','data');