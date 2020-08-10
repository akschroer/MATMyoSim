function demo_twitch_OM
% Function illustrates how to run a simulation of a single-half-sarcomere
% held isometric and activated by a transient pulse of Ca2+

fprintf(1,'\n\n<1 for yes> <0 for no>');
plot_experimental=input('plot experimental values?   ');

% Variables
protocol_file_string = 'protocol_exp_con2.txt';%1s.txt';%
model_parameters_json_file_string = 'twitch_Con_model_exp3.json';%bp
options_file_string = 'twitch_1_options.json';
model_output_file_string = '../../temp\twitch_Con_output.myo';
json_struct_con = loadjson('twitch_Con_model_exp.json');
dist=json_struct_con.MyoSim_model.hs_props.myofilaments.bin_min:json_struct_con.MyoSim_model.hs_props.myofilaments.bin_width:json_struct_con.MyoSim_model.hs_props.myofilaments.bin_max;
pops=size(length(dist),2);

% Make sure the path allows us to find the right files
addpath(genpath('../../../../code'));

% Run a control simulation
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
[mx, ix]=max(vel(500:end));
[mn, in]=min(vel(500:end));
con_time(1) = (ix-in)*ts;
figure(5);
plot(sim_output.time_s,vel,'k-','LineWidth',2)
hold on
plot((500+[ix in])*ts,[mx mn],'g*')
plot((500+[ix in])*ts,[0 0],'LineWidth',2)
figure(6);
hold on
plot(sim_output.time_s,sim_output.f_activated,'k-','LineWidth',2)
plot(sim_output.time_s,sim_output.f_bound,'k--','LineWidth',2)

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output
baseline_force(1) = sim_output.muscle_force(end);
[maf inf]= max(sim_output.muscle_force);%cb_force
peak(1) = maf-baseline_force(1);
AUC(1) = sum(sim_output.muscle_force*ts);
SRX(1,1) = min(sim_output.M1);
SRX(1,2) = sim_output.M1(end);

for t = 1:5
for ib=1:length(dist)
pops(ib) = sim_output.cb_pops(400+200*t,1,ib);
end
figure(7)
subplot(1,2,1)
hold on
plot(dist,pops,'color',[0.9-0.15*t,0.9-0.15*t,0.9-0.15*t],'LineWidth',2);
xlabel('x (nm)')
end
ylabel('bound myosin distribution');
title('Control')
legend('100ms','300ms','500ms','700ms','900ms')




figure(3);
subplot(3,1,1);
plot(sim_output.time_s,sim_output.muscle_force-baseline_force(1),'k-','LineWidth',2);%cb_force
ylabel('Force (N/m^{-2})');
subplot(3,1,2);
plot(sim_output.time_s,sim_output.hs_length,'k-','LineWidth',2);
ylabel('Half-sarcomere length (nm)');
subplot(3,1,3);
plot(sim_output.time_s,sim_output.M1,'k-','LineWidth',2);
ylabel('Myosin in SRX (M_{OFF})');
figure(4)
subplot(1,2,1);
plot(sim_output.time_s,sim_output.M2,'k-','LineWidth',2);
ylabel('Myosin available for binding (M_{ON})');
subplot(1,2,2);
N_unbound = sim_output.f_activated-sim_output.f_bound;
plot(sim_output.time_s,N_unbound,'k-','LineWidth',2);
ylabel('Thin filament available for binding (N_{UNBOUND})');


%%
% Variables_1s.txt';%
protocol_file_string = 'protocol_exp_con2.txt';%_1s';%
model_parameters_json_file_string = 'twitch_OM_model_exp3.json';%_bp
options_file_string = 'twitch_1_options.json';
model_output_file_string = '..\..\temp\twitch_a_output.myo';

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
[mx, ix]=max(vel(500:end));
[mn, in]=min(vel(500:end));
con_time(2) = (ix-in)*ts
baseline_force(2) = sim_output.muscle_force(end);
[maf inf]= max(sim_output.muscle_force);
peak(2) = maf-baseline_force(2)
peak_ratio = peak(2)/peak(1)
AUC(2) = sum(sim_output.muscle_force*ts)
Force_Index = AUC(2)/AUC(1)
SRX(2,1) = min(sim_output.M1);
SRX(2,2) = sim_output.M1(end);
SRX

for t = 1:5
for ib=1:length(dist)
pops(ib) = sim_output.cb_pops(400+200*t,1,ib);
end
figure(7)
subplot(1,2,2)
hold on
plot(dist,pops,'color',[0.9-0.15*t,1,0.9-0.15*t],'LineWidth',2);
xlabel('x (nm)')
end
ylabel('bound myosin distribution');
title('alpha')
legend('100ms','300ms','500ms','700ms','900ms')

figure(5);
hold on
plot(sim_output.time_s,vel,'g-','LineWidth',2)
plot((500+[ix in])*ts,[mx mn],'y*')
plot((500+[ix in])*ts,[50,50],'LineWidth',2)
figure(6);
hold on
plot(sim_output.time_s,sim_output.f_activated,'g-','LineWidth',2)
plot(sim_output.time_s,sim_output.f_bound,'g--','LineWidth',2)
legend('Control activated','Control bound','OM activated','OM bound')
xlabel('Time (ms)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})
ylabel('Actin available')
% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output


figure(3);
subplot(3,1,1);
hold on
plot(sim_output.time_s,sim_output.muscle_force-baseline_force(2),'g-','LineWidth',2);
xlabel('Time (ms)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})
subplot(3,1,2);
hold on
plot(sim_output.time_s,sim_output.hs_length,'g-','LineWidth',2);
xlabel('Time (ms)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})
subplot(3,1,3);
hold on
plot(sim_output.time_s,sim_output.M1,'g-','LineWidth',2);
xlabel('Time (ms)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xticklabels({'-500','0','500','1000','1500'})

figure(4);
subplot(1,2,1);
hold on
plot(sim_output.time_s,sim_output.M2,'g-','LineWidth',2);
subplot(1,2,2);
hold on
N_unbound = sim_output.f_activated-sim_output.f_bound;
plot(sim_output.time_s,N_unbound,'g-','LineWidth',2);




if plot_experimental
    t_exp = load('cell_time_trace.mat');
    Con_F = load('Control_c8_2_force.mat');
    P710R_F = load('P710R_c3_force.mat');
    figure(3)
    subplot(3,1,1);
    hold on
    plot(t_exp.cell_time_trace+0.5,Con_F.Control_c8_2_force,'k:','LineWidth',3);%+113.5
    plot(t_exp.cell_time_trace+0.5,P710R_F.P710R_c3_force,'g:','LineWidth',3);%+1269.1
       legend('Control model','OM model','Control exp.','OM exp.')
    subplot(3,1,3)
    hold on
    plot([0.2 1],[0.43 0.8; 0.43 0.8],'k:')
    plot([0.2 1],[0.1 0.7; 0.1 0.7],'r:')  

end

