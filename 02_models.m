%% This tutorial is largely based on the following CoSMoMVPA tutorial 
% http://www.cosmomvpa.org/matlab/demo_fmri_searchlight_rsm.html#demo-fmri-searchlight-rsm

%% First, let's visualize some models 
% We will use the following data:
% Connolly et al (2012), Representation of biological classes in the
%    human brain. Journal of Neuroscience,
%    doi 10.1523/JNEUROSCI.5547-11.2012
%
%    Six categories (monkey, lemur, mallard, warbler, ladybug, lunamoth)
%    during ten runs in an fMRI study. Using the General Linear Model
%    response were estimated for each category in each run, resulting
%    in 6*10=60 t-values.
ak6_study_path='~\Documents\MATLAB\CoSMoMVPA\tutorial_data\ak6';
output_path='~\Documents\rsa_tutorial\output';

labels={'monkey','lemur',...
        'mallard','warbler',...
        'ladybug','lunamoth'}';
nsamples=length(labels);
animal_class=[1 1 2 2 3 3]';

% 1. Simple RDM 
% similarity structure is perfectly clustered
% primates (monkey, lemur), birds (mallard, warbler)
% and insects (ladybug, lunamoth) are the same (distance=0)
% all other pairs are equally dissimilar (distance=1)

for row=1:nsamples
    for col=1:nsamples
        same_animal_class=animal_class(row)==animal_class(col);
        if same_animal_class
            simple(row,col)=0;
        else
            simple(row,col)=2;
        end
    end
end

figure(1);
rdm_plot1=imagesc(simple,[0,2]);
title('Simple RDM');
set(gca,'XTick',1:nsamples,'XTickLabel',labels,...
        'YTick',1:nsamples,'YTickLabel',labels)
colormap('hot');
colorbar

disp(simple);

% 2. Linear RDM
% similarity structure is linear
% primates of distance 1 from birds and distance 2 from insects
% insects distance 1 from birds
%
% Should use Spearman rather than Pearson correlation measure here

linear=abs(bsxfun(@minus,animal_class,animal_class'));
figure(2);
rdm_plot1=imagesc(linear,[0,2]);
title('Linear RDM');
set(gca,'XTick',1:nsamples,'XTickLabel',labels,...
        'YTick',1:nsamples,'YTickLabel',labels)
colormap('hot');
colorbar

disp(linear);

% 3. Behavioral RDM
% similarity structure is based on 
% behavioral similarity ratings
load(fullfile(ak6_study_path,'models','behav_sim.mat'));

figure(3);
rdm_plot2=imagesc(behav,[0,2]);
title('Behavioral RDM');
set(gca,'XTick',1:nsamples,'XTickLabel',labels,...
        'YTick',1:nsamples,'YTickLabel',labels)
colormap('hot');
colorbar

disp(behav);

%% Now, let's run a single subject RSA searchlight using the RDMs
measure_type=behav;

for i = 1:8
    subjID = strcat('s0',num2str(i));
    
    % set study paths
    data_path=fullfile(ak6_study_path,subjID); % data from subject
    mask_fn=fullfile(data_path, 'brain_mask.nii'); % whole brain mask
    data_fn=fullfile(data_path,'glm_T_stats_even.nii');
    ds=cosmo_fmri_dataset(data_fn,...
                          'mask',mask_fn,...
                          'targets',1:6,...
                          'chunks',1);

    % Set animal species & class
    ds.sa.labels=labels;
    ds.sa.animal_class=animal_class;

    % simple sanity check 
    cosmo_check_dataset(ds);

    % Define Searchlight
    nvoxels_per_searchlight=100; % define neighborhood for each feature (voxel).
    nbrhood=cosmo_spherical_neighborhood(ds,'count',nvoxels_per_searchlight);
    cosmo_disp(nbrhood);

    % set measure
    target_dsm=measure_type;
    measure=@cosmo_target_dsm_corr_measure;
    measure_args=struct();
    measure_args.target_dsm=target_dsm;

    % run searchlight
    ds_rsm_behav=cosmo_searchlight(ds,nbrhood,measure,measure_args);

    % show results
    figure;
    cosmo_plot_slices(ds_rsm_behav);

    % store results
    filename=strcat(subjID,'_even_rsm_linear.nii');
    output_fn=fullfile(output_path,filename);
    cosmo_map2fmri(ds_rsm_behav,output_fn);
end