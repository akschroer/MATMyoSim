function demo_twitch_singles
% Function illustrates how to run a simulation of a single-half-sarcomere
% held isometric and activated by a transient pulse of Ca2+

% Variables
protocol_file_string = 'protocol_exp_con2.txt';%'protocol_1s.txt';%
model_parameters_json_file_string = 'twitch_Con_model.json';
options_file_string = 'twitch_1_options.json';
model_output_file_string = '..\..\temp\twitch_Con_output.myo';
json_struct_con = loadjson('twitch_con_model.json');

% Make sure the path allows us to find the right files
addpath(genpath('..\..\..\..\code'));

% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(1) = (ix-in)*ts;

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output
baselines(1) = sim_output.muscle_force(500);
peak(1) = max(sim_output.muscle_force)-baselines(1);
AUC(1) = sum(sim_output.muscle_force*ts);
SRX(1,1) = min(sim_output.M1);
SRX(2,1) = sim_output.M1(end);
Actin(1) = max(sim_output.f_activated);
%%
figure(3);
subplot(3,1,1)
plot(sim_output.time_s,sim_output.muscle_force-baselines(1),'k-','LineWidth',2);
ylabel('Force (N/m^{-2})');
xlabel('Time (ms)')
axis([0 2 -50 6000])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})

subplot(3,1,2)
plot(sim_output.time_s,sim_output.muscle_force-baselines(1),'k-','LineWidth',2);
ylabel('Force (N/m^{-2})');
xlabel('Time (ms)')
axis([0 2 -50 10000])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})

subplot(3,1,3)
plot(sim_output.time_s,sim_output.muscle_force-baselines(1),'k-','LineWidth',2);
ylabel('Force (N/m^{-2})');
xlabel('Time (ms)')
axis([0 2 -50 10000])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})


%% define additional parameter sets
% Variables
model_output_file_string = '..\..\temp\twitch_P710R_output.myo';
json_struct_710 = loadjson('twitch_P710R_model.json');
json_struct_710_kdet = json_struct_710;
json_struct_710_k3 = json_struct_710;
json_struct_710_k1 = json_struct_710;
json_struct_710_ss = json_struct_710;

json_struct_710_kdetk3 = json_struct_710;
json_struct_710_kdetk1 = json_struct_710;
json_struct_710_k1k3 = json_struct_710;
json_struct_710_k1k3kdet = json_struct_710;

json_struct_710_kdet.MyoSim_model.hs_props.parameters.k_1=json_struct_con.MyoSim_model.hs_props.parameters.k_1;
json_struct_710_kdet.MyoSim_model.hs_props.parameters.k_3=json_struct_con.MyoSim_model.hs_props.parameters.k_3;
json_struct_710_kdet.MyoSim_model.hs_props.parameters.x_ps=json_struct_con.MyoSim_model.hs_props.parameters.x_ps;
jsonStr = jsonencode(json_struct_710_kdet);
fid = fopen('twitch_P710R_model_kdet.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_k3.MyoSim_model.hs_props.parameters.k_1=json_struct_con.MyoSim_model.hs_props.parameters.k_1;
json_struct_710_k3.MyoSim_model.hs_props.parameters.k_4_0=json_struct_con.MyoSim_model.hs_props.parameters.k_4_0;
json_struct_710_k3.MyoSim_model.hs_props.parameters.k_4_1=json_struct_con.MyoSim_model.hs_props.parameters.k_4_1;
json_struct_710_k3.MyoSim_model.hs_props.parameters.x_ps=json_struct_con.MyoSim_model.hs_props.parameters.x_ps;
jsonStr = jsonencode(json_struct_710_k3);
fid = fopen('twitch_P710R_model_k3.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_k1.MyoSim_model.hs_props.parameters.k_3=json_struct_con.MyoSim_model.hs_props.parameters.k_3;
json_struct_710_k1.MyoSim_model.hs_props.parameters.k_4_0=json_struct_con.MyoSim_model.hs_props.parameters.k_4_0;
json_struct_710_k1.MyoSim_model.hs_props.parameters.k_4_1=json_struct_con.MyoSim_model.hs_props.parameters.k_4_1;
json_struct_710_k1.MyoSim_model.hs_props.parameters.x_ps=json_struct_con.MyoSim_model.hs_props.parameters.x_ps;
jsonStr = jsonencode(json_struct_710_k1);
fid = fopen('twitch_P710R_model_k1.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_ss.MyoSim_model.hs_props.parameters.k_3=json_struct_con.MyoSim_model.hs_props.parameters.k_3;
json_struct_710_ss.MyoSim_model.hs_props.parameters.k_4_0=json_struct_con.MyoSim_model.hs_props.parameters.k_4_0;
json_struct_710_ss.MyoSim_model.hs_props.parameters.k_4_1=json_struct_con.MyoSim_model.hs_props.parameters.k_4_1;
json_struct_710_ss.MyoSim_model.hs_props.parameters.k_1=json_struct_con.MyoSim_model.hs_props.parameters.k_1;
jsonStr = jsonencode(json_struct_710_ss);
fid = fopen('twitch_P710R_model_ss.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);

json_struct_710_kdetk3.MyoSim_model.hs_props.parameters.k_1=json_struct_con.MyoSim_model.hs_props.parameters.k_1;
jsonStr = jsonencode(json_struct_710_kdetk3);
fid = fopen('twitch_P710R_model_kdetk3.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_kdetk3.MyoSim_model.hs_props.parameters.k_3=json_struct_con.MyoSim_model.hs_props.parameters.k_3;
jsonStr = jsonencode(json_struct_710_kdetk1);
fid = fopen('twitch_P710R_model_kdetk1.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_k1k3.MyoSim_model.hs_props.parameters.k_4_0=json_struct_con.MyoSim_model.hs_props.parameters.k_4_0;
json_struct_710_k1k3.MyoSim_model.hs_props.parameters.k_4_1=json_struct_con.MyoSim_model.hs_props.parameters.k_4_1;
jsonStr = jsonencode(json_struct_710_k1k3);
fid = fopen('twitch_P710R_model_k1k3.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
json_struct_710_k1k3kdet.MyoSim_model.hs_props.parameters.x_ps=json_struct_con.MyoSim_model.hs_props.parameters.x_ps;
jsonStr = jsonencode(json_struct_710_k1k3kdet);
fid = fopen('twitch_P710R_model_k1k3kdet.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);

%% run with only k1

model_parameters_json_file_string = 'twitch_P710R_model_k1.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(2) = (ix-in)*ts;
baselines(2) = sim_output.muscle_force(500);
peak(2) = max(sim_output.muscle_force)-baselines(2);
AUC(2) = sum(sim_output.muscle_force*ts);
SRX(1,2) = min(sim_output.M1);
SRX(2,2) = sim_output.M1(end);
Actin(2) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,2)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(2),'r--','LineWidth',1);

%% run with only k3

model_parameters_json_file_string = 'twitch_P710R_model_k3.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(3) = (ix-in)*ts;
baselines(3) = sim_output.muscle_force(500);
peak(3) = max(sim_output.muscle_force)-baselines(3);
AUC(3) = sum(sim_output.muscle_force*ts);
SRX(1,3) = min(sim_output.M1);
SRX(2,3) = sim_output.M1(end);
Actin(3) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,2)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(3),'r-.','LineWidth',1);

%% run with only step size

model_parameters_json_file_string = 'twitch_P710R_model_ss.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(4) = (ix-in)*ts;
baselines(4) = sim_output.muscle_force(500);
peak(4) = max(sim_output.muscle_force)-baselines(4);
AUC(4) = sum(sim_output.muscle_force*ts);
SRX(1,4) = min(sim_output.M1);
SRX(2,4) = sim_output.M1(end);
Actin(4) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,2)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(4),'r:','LineWidth',2);

%% run with only kdet

model_parameters_json_file_string = 'twitch_P710R_model_kdet.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(5) = (ix-in)*ts;
baselines(5) = sim_output.muscle_force(500);
peak(5) = max(sim_output.muscle_force)-baselines(5);
AUC(5) = sum(sim_output.muscle_force*ts);
SRX(1,5) = min(sim_output.M1);
SRX(2,5) = sim_output.M1(end);
Actin(5) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,2)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(5),'r:','LineWidth',1);

legend('control','only k1','only k3','only stepsize','only kdet')
%[0 1.2 0 250000]

%% run missing k1

model_parameters_json_file_string = 'twitch_P710R_model_kdetk3.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(6) = (ix-in)*ts;
baselines(6) = sim_output.muscle_force(500);
peak(6) = max(sim_output.muscle_force)-baselines(6);
AUC(6) = sum(sim_output.muscle_force*ts);
SRX(1,6) = min(sim_output.M1);
SRX(2,6) = sim_output.M1(end);
Actin(6) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,3)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(6),'r--','LineWidth',2);


%% run missing k3

model_parameters_json_file_string = 'twitch_P710R_model_kdetk1.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(7) = (ix-in)*ts;
baselines(7) = sim_output.muscle_force(500);
peak(7) = max(sim_output.muscle_force)-baselines(7);
AUC(7) = sum(sim_output.muscle_force*ts);
SRX(1,7) = min(sim_output.M1);
SRX(2,7) = sim_output.M1(end);
Actin(7) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,3)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(7),'r-.','LineWidth',2);



%% run missing stepsize

model_parameters_json_file_string = 'twitch_P710R_model_k1k3kdet.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(8) = (ix-in)*ts;
baselines(8) = sim_output.muscle_force(500);
peak(8) = max(sim_output.muscle_force)-baselines(8);
AUC(8) = sum(sim_output.muscle_force*ts);
SRX(1,8) = min(sim_output.M1);
SRX(2,8) = sim_output.M1(end);
Actin(8) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,3)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(8),'r:','LineWidth',3);
%[0 1.2 0 250000]

%% run missing kdet

model_parameters_json_file_string = 'twitch_P710R_model_k1k3.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(9) = (ix-in)*ts;
baselines(9) = sim_output.muscle_force(500);
peak(9) = max(sim_output.muscle_force)-baselines(9);
AUC(9) = sum(sim_output.muscle_force*ts);
SRX(1,9) = min(sim_output.M1);
SRX(2,9) = sim_output.M1(end);
Actin(9) = max(sim_output.f_activated);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output
%%
figure(3);
subplot(3,1,3)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(9),'r:','LineWidth',2);
legend('control','not k1','not k3','not stepsize','not kdet')
axis([0 2 0 10000])

%% run with all

model_parameters_json_file_string = 'twitch_P710R_model.json';
% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_json_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

len = length(sim_output.hs_length);
vel = zeros(len,1);
ts = sim_output.time_s(2)-sim_output.time_s(1);
for i=2:len
vel(i)=(sim_output.hs_length(i)-sim_output.hs_length(i-1))/ts;
end
[mx, ix]=max(vel(400:end));
[mn, in]=min(vel(400:end));
con_time(10) = (ix-in)*ts
baselines(10) = sim_output.muscle_force(500);
peak(10) = max(sim_output.muscle_force)-baselines(10)
peak_ratios = peak/peak(1)
AUC(10) = sum(sim_output.muscle_force*ts)
Force_Index = AUC/AUC(1)
SRX(1,10) = min(sim_output.M1);
SRX(2,10) = sim_output.M1(end)
Actin(10) = max(sim_output.f_activated)

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(3);
subplot(3,1,1)
hold on
plot(sim_output.time_s,sim_output.muscle_force-baselines(10),'r-','LineWidth',2);
legend('control','P710R')
axis([0 2 0 6000])%[0 1.2 0 250000]


