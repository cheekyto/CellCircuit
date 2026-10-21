%% Load data
n_replicates = size(dato,1);
netBurst_peakFiringNorm = nan(2,7,n_replicates);

% baseline
netBurst_peakFiringNorm(1,1,:) = dato(:,1);   % WT
netBurst_peakFiringNorm(2,1,:) = dato(:,2);   % het
% gabazine 10 uM
netBurst_peakFiringNorm(1,2,:) = dato(:,3);
netBurst_peakFiringNorm(2,2,:) = dato(:,4);
% gabazine 20 uM
netBurst_peakFiringNorm(1,3,:) = dato(:,5);
netBurst_peakFiringNorm(2,3,:) = dato(:,6);
% CGP 10 uM
netBurst_peakFiringNorm(1,4,:) = dato(:,7);
netBurst_peakFiringNorm(2,4,:) = dato(:,8);
% CGP 20 uM
netBurst_peakFiringNorm(1,5,:) = dato(:,9);
netBurst_peakFiringNorm(2,5,:) = dato(:,10);
% D-AP5 
netBurst_peakFiringNorm(1,6,:) = dato(:,11);
netBurst_peakFiringNorm(2,6,:) = dato(:,12);
% DNQX
netBurst_peakFiringNorm(1,7,:) = dato(:,13);
netBurst_peakFiringNorm(2,7,:) = dato(:,14);

switch n_replicates
    case 10
         fprintf("\nData from the project, 105 div \n\n")
    case 13
         fprintf("\nData from the project, 166 div \n\n")
    case 14
         fprintf("\nData from the project, 240 div \n\n")
end

%% Save data

main_folder = '/path/to/main_folder';
project_folder = 'project';
switch n_replicates
    case 10
         group = 'project_105div';
    case 13
         group = 'project_166div';
    case 14
         group = 'project_240div';
end
file_output = sprintf('%s\\%s\\%s\\metrics_%s.mat', main_folder, project_folder, group,group);

save(file_output, 'netBurst_peakFiringNorm','-append')
% spikeRate, entropy, burst_duration, burst_amplitude,
% netBurst_correlation, netBurst_peakSTD, netBurst_onsetSTD,
% netBurst_distanceCV, netBurst_distance, netBurst_participation,
% netBurst_neurons, netBurst_spikes, netBurst_partCV, netBurst_rate,
% netBurst_rateCV, netBurst_duration, netBurst_amplitude,
% netBurst_peakFiring, netBurst_peakFiringNorm

%% Normalize data

metric = spikeRate;

n_replicates = size(metric,3);
metric_norm = nan(2,7,n_replicates);

% baseline
metric_norm(1,1,:) = metric(1,1,:)./metric(1,1,:);   % WT
metric_norm(2,1,:) = metric(2,1,:)./metric(2,1,:);   % het
% gabazine 10 uM
metric_norm(1,2,:) = metric(1,2,:)./metric(1,1,:);
metric_norm(2,2,:) = metric(2,2,:)./metric(2,1,:);
% gabazine 20 uM
metric_norm(1,3,:) = metric(1,3,:)./metric(1,1,:);
metric_norm(2,3,:) = metric(2,3,:)./metric(2,1,:);
% CGP 10 uM
metric_norm(1,4,:) = metric(1,4,:)./metric(1,1,:);
metric_norm(2,4,:) = metric(2,4,:)./metric(2,1,:);
% CGP 20 uM
metric_norm(1,5,:) = metric(1,5,:)./metric(1,1,:);
metric_norm(2,5,:) = metric(2,5,:)./metric(2,1,:);
% D-AP5 
metric_norm(1,6,:) = metric(1,6,:)./metric(1,1,:);
metric_norm(2,6,:) = metric(2,6,:)./metric(2,1,:);
% DNQX
metric_norm(1,7,:) = metric(1,7,:)./metric(1,1,:);
metric_norm(2,7,:) = metric(2,7,:)./metric(2,1,:);

datos_wt = nan(size(metric,3),6);
datos_het = nan(size(metric,3),6);
for j=1:6
    datos_wt(:,j) =  metric_norm(1,j+1,:);
    datos_het(:,j) =  metric_norm(2,j+1,:);
end

%% Plot data

n_replicates = size(metric,3);
figure

for k=1:2
    switch k
        case 1
            label_group = "wt";
        case 2
            label_group = "het";
    end

    subplot(1,2,k)

    label = "baseline";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,1,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y
    ylim([0 max(max(max(metric)))])

    label = "gabazine 10";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,2,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "gabazine 20";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,3,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "CGP 10";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,4,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "CGP 20";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,5,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

if n_replicates==13
    label = "DNQX";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,7,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y
end
    title(label_group)
end

metric = metric_norm;
figure
for k=1:2
    switch k
        case 1
            label_group = "wt";
        case 2
            label_group = "het";
    end

    subplot(1,2,k)

    label = "baseline";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,1,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y
    ylim([0 max(max(max(metric)))])

    label = "gabazine 10";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,2,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "gabazine 20";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,3,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "CGP 10";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,4,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

    label = "CGP 20";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,5,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y

if n_replicates==13
    label = "DNQX";
    for i=1:n_replicates
        labelx(i) = label;
    end
    y = nan(n_replicates,1);
    y(:) = metric(k,7,:);
    swarmchart(categorical(labelx), y,'XJitterWidth', 0.2)
    hold on
    bar(categorical(label), median(y,'omitnan'), 'FaceAlpha',0.5)
    clear labelx y
end
    title(label_group)
end


%% 

% spikeRate, entropy, burst_duration, burst_amplitude,
% netBurst_correlation, netBurst_peakSTD, netBurst_onsetSTD,
% netBurst_distanceCV, netBurst_distance, netBurst_participation,
% netBurst_neurons, netBurst_spikes, netBurst_partCV, netBurst_rate,
% netBurst_rateCV, netBurst_duration, netBurst_amplitude
% netBurst_peakFiring, netBurst_peakFiringNorm

metric1 = netBurst_correlation;
metric2 = netBurst_peakFiringNorm;
n_cond = 5;

n_replicates = size(metric1,3);
aux1 = nan(2*n_cond,1);
aux2 = nan(2*n_cond,1);

idx=0;
for i=1:2
    for j=1:n_cond
        for k=1: n_replicates
            idx=idx+1;
            aux1(idx,1) = metric1(i,j,k);
            aux2(idx,1) = metric2(i,j,k);
        end
    end
end

figure 
plot(aux1, aux2, '*')