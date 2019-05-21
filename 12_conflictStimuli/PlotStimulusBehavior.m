% Fish{1} = '/Users/migault/PhD/Pr?sentations/2019/JournalClub/Behavior2019-05-15/F1';
% Fish{2} = '/Users/migault/PhD/Pr?sentations/2019/JournalClub/Behavior2019-05-15/F2';
% 
% FishTitle{1} = 'Fish 1 : 2019/05/15';
% FishTitle{1} = 'Fish 1 : 2019/05/15';
% 
% for i = 1:size(Fish, 2)
%     for j = [5, 10, 20]
%         Mot = readtable([Fish{i} '/', num2str(j), '/Stimulus.txt'],'Format','%f%f%f');
%         Vis = readtable([Fish{i} '/', num2str(j), '/StimulusVisualVestibular.txt'],'Format','%f%f%f');
%         subplot(size(Fish, 2), j, j+(i-1)*3)
%         plot(Mot.Var2-Mot.Var2(1), Mot.Var3);
%         hold on
%         plot(Vis.Time-Vis.Time(1), Vis.Visual);
%     end
% end

% Define path to run and to tracking file (DO NOT FORGET TO CHANGE THAT):
run_folder = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-05-15_Behavior/Run 15';
tracking = '15h46m21s-Tracking.txt';

% Get info from .txt files:
Mot = readtable(fullfile(run_folder, 'Stimulus.txt'), 'Format', '%f%f%f');
Vis = readtable(fullfile(run_folder, 'StimulusVisualVestibular.txt'), 'Format', '%f%f%f');
Tra = readtable(fullfile(run_folder, tracking), 'Format', '%f%f%f');

% Split into variables:
time1 = Vis.Time-Vis.Time(1);
motor = Vis.Motor;
visual = Vis.Visual ./ max(Vis.Visual) * max(motor);
time2 = Tra.Var1-Tra.Var1(1);
tracking_brut = Tra.Var2;
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.15,'DesignMethod','butter');
tracking = filtfilt(d1, tracking_brut);
tracking_nobl = tracking - baselinefit(tracking')';
track_interp = interp1(time2, tracking_nobl, time1);

% Plotting:
figure
subplot(3, 1, 1)
hold on
plot(time1, motor);
plot(time1, visual);
subplot(3, 1, 2)
hold on
plot(time2, tracking_brut, '-.')
plot(time2, tracking)
subplot(3, 1, 3)
plot(time2, tracking_nobl, 'k')

% Defining variables:
var1 = motor;
var2 = visual;
var3 = ~isnan(var2);
var2(isnan(var2)) = 1;
var4 = var1 + var2;

% Multilinear regression:
vartot = [ones(size(var1)), var1, var2, var3, var4];
[b, ~, r, ~, stats] = regress(track_interp, vartot);
figure; 
hold on
plot(time1, track_interp)
plot(time1, sum(b'.*vartot, 2))





