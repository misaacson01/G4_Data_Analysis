function G4_Plot_Data_flyingdetector(exp_folder, trial_options, CL_conds, OL_conds, TC_conds)
%FUNCTION G4_Plot_Data_flyingdetector(exp_folder, trial_options, CL_conds, OL_conds, TC_conds)
% 
% Inputs:
% exp_folder: path containing G4_Processed_Data.mat file
% trial_options: 1x3 logical array [pre-trial, intertrial, post-trial]
% CL_conds: matrix of closed-loop (CL) conditions to plot as histograms
% OL_conds: matrix of open-loop (OL) conditions to plot as timeseries
% TC_conds: matrix of open-loop conditions to plot as tuning curves (TC)


%% user-defined parameters
%datatype options: 'LmR_chan', 'L_chan', 'R_chan', 'F_chan', 'Frame Position', 'LmR', 'LpR'
CL_datatypes = {'Frame Position'}; %datatypes to plot as histograms
OL_datatypes = {'LmR','LpR'}; %datatypes to plot as timeseries
TC_datatypes = {'LmR','LpR'}; %datatypes to plot as tuning curves

%specify plot properties
rep_Color = [0.5 0.5 0.5];
mean_Color = [0 0 0];
rep_LineWidth = 0.05;
mean_LineWidth = 1;
subtitle_FontSize = 8;


%% load data and prepare for plotting
%load G4_Processed_Data
if nargin==0
    exp_folder = uigetdir('C:/','Select a folder containing a G4_Processed_Data file');
    trial_options = [0 0 0];
end
files = dir(exp_folder);
try
    Data_name = files(contains({files.name},{'G4_Processed_Data'})).name;
catch
    error('cannot find G4_Processed_Data file in specified folder')
end
load(fullfile(exp_folder,Data_name));

%create default matrices for plotting all conditions
if nargin<5 
    CL_conds = find(Data.conditionModes==4); %all closed-loop modes
    OL_conds = find(Data.conditionModes~=4); %all open-loop modes
    TC_conds = []; %by default, don't plot any tuning curves
    CL_conds = reshape(CL_conds,[round(sqrt(numel(CL_conds))) ceil(sqrt(numel(CL_conds)))])';
    OL_conds = reshape(OL_conds,[round(sqrt(numel(OL_conds))) ceil(sqrt(numel(OL_conds)))])';
end

%get datatype indices
num_OL_datatypes = length(OL_datatypes);
OL_inds = nan(1,num_OL_datatypes);
for i = 1:num_OL_datatypes
    OL_inds(i) = find(strcmpi(Data.channelNames.timeseries,OL_datatypes{i}));
end
num_CL_datatypes = length(CL_datatypes);
CL_inds = nan(1,num_CL_datatypes);
for i = 1:num_CL_datatypes
    CL_inds(i) = find(strcmpi(Data.channelNames.timeseries,CL_datatypes{i}));
end
num_TC_datatypes = length(TC_datatypes);
TC_inds = nan(1,num_TC_datatypes);
for i = 1:num_TC_datatypes
    TC_inds(i) = find(strcmpi(Data.channelNames.timeseries,TC_datatypes{i}));
end


%% plot data
%calculate overall measurements and plot basic histograms
figure()
for i = 1:num_TC_datatypes
    subplot(2,num_TC_datatypes,1)
    data_vec = reshape(Data.timeseries(TC_inds(i),:,:,:),[1 numel(Data.timeseries(TC_inds(i),:,:,:))]);
    text(0.1, 0.9-0.2*i, ['Mean ' TC_datatypes{i} ' = ' num2str(nanmean(data_vec))]);
    axis off
    hold on
    
    subplot(2,num_TC_datatypes,num_TC_datatypes+i)
    avg = length(data_vec)/100;
    hist(data_vec,100)
    hold on
    xl = xlim;
    plot(xl,[avg avg],'--','Color',rep_Color','LineWidth',mean_LineWidth)
    title([TC_datatypes{i} ' Histogram'],'FontSize',subtitle_FontSize);
end
if trial_options(2)==1
    subplot(2,num_TC_datatypes,num_TC_datatypes)
    plot(Data.interhistogram','Color',rep_Color,'LineWidth',rep_LineWidth)
    hold on
    plot(nanmean(Data.interhistogram),'Color',mean_Color,'LineWidth',mean_LineWidth)
    title('Intertrial Pattern Frame Histogram','FontSize',subtitle_FontSize)
end

%plot histograms for closed-loop trials
if ~isempty(CL_conds)
    num_figs = size(CL_conds,3);
    for d = 1:length(CL_datatypes)
        for fig = 1:num_figs
            num_plot_rows = max(nansum(CL_conds(:,:,fig)>0));
            num_plot_cols = max(nansum(CL_conds(:,:,fig)>0,2));
            figure()
            for row = 1:num_plot_rows
                for col = 1:num_plot_cols
                    cond = CL_conds(row,col,fig);
                    if cond>0
                        subplot(num_plot_rows, num_plot_cols, col+num_plot_cols*(row-1))
                        [~, ~, num_reps, num_positions] = size(Data.histograms);
                        x = circshift(1:num_positions,[1 num_positions/2]);
                        x(x>x(end)) = x(x>x(end))-num_positions;
                        tmpdata = circshift(squeeze(Data.histograms(d,cond,:,:)),[1 num_positions/2]);
                        plot(repmat(x',[1 num_reps]),tmpdata','Color',rep_Color,'LineWidth',rep_LineWidth);
                        hold on
                        plot(x,nanmean(tmpdata),'Color',mean_Color,'LineWidth',mean_LineWidth)
                        title(['Condition #' num2str(cond)],'FontSize',subtitle_FontSize)
                    end
                end
            end
        end
    end
end

%plot timeseries data for open-loop trials
if ~isempty(OL_conds)
    num_figs = size(OL_conds,3);
    num_reps = size(data.timeseries,3);
    %loop for different data types
    for d = 1:length(OL_datatypes)
        for fig = 1:num_figs
            num_plot_rows = max(nansum(OL_conds(:,:,fig)>0));
            num_plot_cols = max(nansum(OL_conds(:,:,fig)>0,2));
            figure()
            for row = 1:num_plot_rows
                for col = 1:num_plot_cols
                    subplot(num_plot_rows, num_plot_cols, col+num_plot_cols*(row-1))
                    subplot(num_plot_rows, num_plot_cols, col+num_plot_cols*(row-1))
                    plot(repmat(Data.timestamps',[1 num_reps]),squeeze(Data.timeseries(d,cond,:,:))','Color',rep_Color,'LineWidth',rep_LineWidth);
                    hold on
                    plot(Data.timestamps,squeeze(nanmean(Data.timeseries(s,cond,:,:),3)),'Color',mean_Color,'LineWidth',mean_LineWidth);
                    title(['Condition #' num2str(cond)],'FontSize',subtitle_FontSize)
                end
            end
        end
    end
end

%plot tuning-curves for specified open-loop trials
if ~isempty(TC_conds)
    num_figs = size(TC_conds,3);
    %loop for different data types
    for d = 1:length(TC_datatypes)
        for fig = 1:num_figs
            num_plot_rows = max(nansum(TC_conds(:,:,fig)>0));
            figure()
            for row = 1:num_plot_rows
                conds = TC_conds(row,:,fig);
                conds(isnan(conds)) = [];
                subplot(num_plot_rows, 1, row)
                plot(squeeze(Data.summaries(d,conds,:)),'Color',rep_Color,'LineWidth',rep_LineWidth);
                hold on
                plot(nanmean(Data.summaries(d,conds,:),3),'Color',mean_Color,'LineWidth',mean_LineWidth);
                title(['Condition #' num2str(cond)],'FontSize',subtitle_FontSize)
            end
        end
    end
end