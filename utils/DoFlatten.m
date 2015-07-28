function[] = DoFlatten(data, filename, type, res)

switch res
    case '32'
        a = gifti('test4.gii');
        a.cdata = squeeze(data');
        save(a,[filename '.gii'],'ExternalFileBinary');
    case '164'
        sub = 100307;
        subDir = ['/a/documents/connectome/_all/' num2str(sub) '/'];
        a = gifti([subDir 'MNINonLinear/' num2str(sub) '.L.SmoothedMyelinMap.164k_fs_LR.func.gii']);
        a.cdata = squeeze(data');
        save(a,[filename '.gii'],'ExternalFileBinary');
end

switch type
    case 'label'
        % create label txt
        num=length(nonzeros(unique(data)));
        cols = cbrewer('qual', 'Set1', num);
        fid = fopen('label_text.txt','w');
        for i = 1:num
            fprintf(fid, '%s \n%3i %3i %3i %3i %3i\n', ...
                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
        end
        fclose(fid);

        unix(['wb_command -metric-label-import ' filename '.gii label_text.txt ' filename '.label']);
        unix(['wb_command -cifti-create-label ' filename '.dlabel.nii -left-label ' filename '.label']);

    case 'label_col_left'
        % create label txt
        %num=length(nonzeros(unique(data)));
        num=17;
        %cols = cbrewer('qual', 'Set1', num);
        t = [0 139 159; 255 146 0; 47 165 80 ]./255;
R = [linspace(t(1,1),t(2,1),9) linspace(t(2,1),t(3,1),9)];
G = [linspace(t(1,2),t(2,2),9) linspace(t(2,2),t(3,2),9)];
B = [linspace(t(1,3),t(2,3),9) linspace(t(2,3),t(3,3),9)];
cols = [R([1:9 11:18])', G([1:9 11:18])', B([1:9 11:18])'];
%cols = [0 0 0; cols];
        %cols = [cbrewer('qual', 'Set1', num - 1); 1 1 1 ];
%        cols = ([0 139 159; 255 146 0; 47 165 80])./255;
        fid = fopen('label_text.txt','w');
        for i = 1:num
%            fprintf(fid, '%s \n%3i %3i %3i %3i %3i\n', ...
%                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
            fprintf(fid, '%s \n%3i %3.0f %3.0f %3.0f %3i\n', ...
                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
        end
        fclose(fid);

        unix(['wb_command -metric-label-import ' filename '.gii label_text.txt ' filename '.label']);
        unix(['wb_command -cifti-create-label ' filename '.dlabel.nii -left-label ' filename '.label']);
    case 'label_col_right'
        % create label txt
%         num=length(nonzeros(unique(data)));
%         cols = cbrewer('qual', 'Set1', num);
%         %cols = cbrewer('seq', 'YlOrRd', num);
%         % cols = [cbrewer('qual', 'Set1', num - 1); 1 1 1 ];
% %        cols = ([0 139 159; 255 146 0; 47 165 80])./255;
%         fid = fopen('label_text.txt','w');
%         for i = 1:num
%             fprintf(fid, '%s \n%3i %3i %3i %3i %3i\n', ...
%                 ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
%         end
%         fclose(fid);

        num=17;
        %cols = cbrewer('qual', 'Set1', num);
        t = [0 139 159; 255 146 0; 47 165 80 ]./255;
R = [linspace(t(1,1),t(2,1),9) linspace(t(2,1),t(3,1),9)];
G = [linspace(t(1,2),t(2,2),9) linspace(t(2,2),t(3,2),9)];
B = [linspace(t(1,3),t(2,3),9) linspace(t(2,3),t(3,3),9)];
cols = [R([1:9 11:18])', G([1:9 11:18])', B([1:9 11:18])'];
%cols = [0 0 0; cols];
        %cols = [cbrewer('qual', 'Set1', num - 1); 1 1 1 ];
%        cols = ([0 139 159; 255 146 0; 47 165 80])./255;
        fid = fopen('label_text.txt','w');
        for i = 1:num
%            fprintf(fid, '%s \n%3i %3i %3i %3i %3i\n', ...
%                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
            fprintf(fid, '%s \n%3i %3.0f %3.0f %3.0f %3i\n', ...
                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
        end
        fclose(fid);
        

        unix(['wb_command -metric-label-import ' filename '.gii label_text.txt ' filename '.label']);
        unix(['wb_command -cifti-create-label ' filename '.dlabel.nii -right-label ' filename '.label']);
        
    case 'scalar_left'
        unix(['wb_command -cifti-create-dense-scalar ' filename '.dscalar.nii -left-metric ' filename '.gii']);

    case 'scalar_right'
        unix(['wb_command -cifti-create-dense-scalar ' filename '.dscalar.nii -right-metric ' filename '.gii']);
       
    case 'label164'
        sub = 100307;
        subDir = ['/a/documents/connectome/_all/' num2str(sub) '/'];
        a = gifti([subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.L.BA.32k_fs_LR.label.gii']);
        a.cdata = squeeze(data');
        save(a,[filename '.gii'],'ExternalFileBinary');
         % create label txt
        num=length(nonzeros(unique(data)));
        cols = cbrewer('qual', 'Dark2', num);
        fid = fopen('label_text.txt','w');
        for i = 1:num
            fprintf(fid, '%s \n%3i %3i %3i %3i %3i\n', ...
                ['Network_' num2str(i)], i, 255 .* cols(i,:), 255);
        end
        fclose(fid);

        unix(['wb_command -metric-label-import ' filename '.gii label_text.txt ' filename '.label']);
        unix(['wb_command -cifti-create-label ' filename '.dlabel.nii -left-label ' filename '.label ' ...
            '-right-label ' subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.R.BA.32k_fs_LR.label.gii']);
        
        unix(['/scr/litauen1/connectome_wb/workbench/bin_linux64/wb_command -cifti-resample ' ...
            filename '.dlabel.nii COLUMN ' ...
            subDir 'MNINonLinear/' num2str(sub) '.BA.164k_fs_LR.dlabel.nii COLUMN ' ...
            'ADAP_BARY_AREA ENCLOSING_VOXEL ' ...
            filename '.164k.dlabel.nii '...
            '-left-spheres ' subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
            '.L.sphere.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.L.sphere.164k_fs_LR.surf.gii '...
            '-left-area-surfs ' ...
            subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.L.midthickness.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.L.midthickness.164k_fs_LR.surf.gii ' ...
            '-right-spheres ' subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
            '.R.sphere.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.R.sphere.164k_fs_LR.surf.gii '...
            '-right-area-surfs ' ...
            subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.R.midthickness.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.R.midthickness.164k_fs_LR.surf.gii ']);
        
    case 'scalar164'
        sub = 100307;
        subDir = ['/a/documents/connectome/_all/' num2str(sub) '/'];
        a = gifti([subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.L.MyelinMap.32k_fs_LR.func.gii']);
        a.cdata = squeeze(data');
        save(a,[filename '.gii'],'ExternalFileBinary');
        
        a = gifti([subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.R.MyelinMap.32k_fs_LR.func.gii']);
        a.cdata = squeeze(data');
        save(a,[filename '.R.gii'],'ExternalFileBinary');
        
        unix(['wb_command -cifti-create-dense-scalar ' filename '.dscalar.nii -left-metric ' filename '.gii']);
        
        unix(['/scr/litauen1/connectome_wb/workbench/bin_linux64/wb_command -cifti-resample ' ...
            filename '.dscalar.nii COLUMN ' ...
            subDir 'MNINonLinear/' num2str(sub) '.MyelinMap.164k_fs_LR.dscalar.nii COLUMN ' ...
            'ADAP_BARY_AREA ENCLOSING_VOXEL ' ...
            filename '.164k.dscalar.nii '...
            '-left-spheres ' subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
            '.L.sphere.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.L.sphere.164k_fs_LR.surf.gii '...
            '-left-area-surfs ' ...
            subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.L.midthickness.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.L.midthickness.164k_fs_LR.surf.gii ' ...
            '-right-spheres ' subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
            '.R.sphere.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.R.sphere.164k_fs_LR.surf.gii '...
            '-right-area-surfs ' ...
            subDir 'MNINonLinear/fsaverage_LR32k/' num2str(sub) '.R.midthickness.32k_fs_LR.surf.gii ' ...
            subDir 'MNINonLinear/' num2str(sub) '.R.midthickness.164k_fs_LR.surf.gii ']);
        

end

% addpath(path, 'toolbox_graph/toolbox/');
% addpath(path, 'toolbox_graph/off/');
% clear options;
% 
% figure;
% vertex=surfm.coord;
% faces=double(surfm.tri);
% plot_mesh(vertex, faces);
% shading interp;
% 
% % you can try with other boundary type
% options.boundary = 'circle';
% % you can try with other Laplacians
% options.laplacian = 'conformal';
% 
% % compute the parameterization
% options.method = 'freeboundary';
% xy1 = compute_parameterization(vertex,faces,options);
% % display the parameterization
% clf;
% faces = double(faces);
% A = triangulation2adjacency(faces);
% plot_graph(A,xy1,'b-'); axis tight;
% 
% options.face_vertex_color = 


