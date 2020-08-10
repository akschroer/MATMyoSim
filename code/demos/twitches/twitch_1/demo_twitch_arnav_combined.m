function demo_twitch_2
% Function illustrates how to run a simulation of a single-half-sarcomere
% held isometric and activated by a transient pulse of Ca2+

%fprintf(1,'\n\n<1 for yes> <0 for no>');
%plot_experimental=input('plot experimental values?   ');
%forced it to be 1 so it does not prompt 
plot_experimental = 1;


outputdir = "./figuresanalysis/"; % path where the output directories will go

% input the file with a list of fully qualified paramters
parameter_file = input('Enter parameter filename: ', 's');

parfid = fopen(parameter_file);
pline = fgetl(parfid);

key_indx = 0;
% loop through all the parameters
while ischar(pline)
    disp(pline)
    


% input the parameter to modify in the model json file and the factor to
% change it by
%create a directory to store the grpahs for the parameter

%parameter_str=input('Enter parameter name: ', 's');

parameter_str = pline;

% track the indx of the parameter    
key_indx = key_indx+1;

% store the name of the parameter 
keylist{key_indx} = strrep(parameter_str,' ','.');

param_indx = 0;
% vary the parameter by the 4 factors 
for parameter_factor = [0.1 0.5 1.0 2.0 10.0]
param_indx = param_indx+1;
dirname = strcat(outputdir,strrep(parameter_str,' ','.'),'-',num2str(parameter_factor),'/');
fprintf(1,"%s \n",dirname);
mkdir(dirname);

% Make sure the path allows us to find the right files
addpath(genpath('../utilities/'));
addpath(genpath('../'));
addpath(genpath('../utilities/jsonlab_1_5/'));
addpath(genpath('../utilities/legendflex/'));

% Variables
protocol_file_string = 'protocol_exp_con2.txt';%1s.txt';%
%model_parameters_json_file_string = 'twitch_Con_model_exp.json';%exp
options_file_string = 'twitch_1_options.json';
model_output_file_string = './temp/twitch_Con_output.myo';
json_struct_con = loadjson('twitch_1_model.json');
dist=json_struct_con.MyoSim_model.hs_props.myofilaments.bin_min:json_struct_con.MyoSim_model.hs_props.myofilaments.bin_width:json_struct_con.MyoSim_model.hs_props.myofilaments.bin_max;
pops=size(length(dist),2);

%change the json parameters
parameter_var=matlab.lang.makeValidName(split(parameter_str));
json_struct_con_exp = loadjson('twitch_Con_model_exp.json');%exp
if(length(parameter_var) == 2)
json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}) = (json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}))*parameter_factor;
end
if(length(parameter_var) == 3)
json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}).(parameter_var{3}) = (json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}).(parameter_var{3}))*parameter_factor;
end

savejson('',json_struct_con_exp,'FileName','twitch_Con_model_exp_changed2.json')
model_parameters_json_file_string = 'twitch_Con_model_exp_changed2.json';%exp

% 
%jsonText = fileread('twitch_Con_model_exp.json');

%jsonData = jsondecode(jsonText); 
% Change HighPrice value in Row 3 from 10000 to 12000
%jsonData{3}.HighPrice = 12000;
% Convert to JSON text
%jsonText2 = jsonencode(jsonData);
% Write to a json file
%fid = fopen('Portfolio2.json', 'w');
%fprintf(fid, '%s', jsonText2);
%fclose(fid);

%%%%%


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
[mx, ix]=max(vel);
[mn, in]=min(vel);
con_time(1) = (ix-in)*ts;

% commented to test without graphs
figure(5);
plot(sim_output.time_s,vel,'k-','LineWidth',2)
hold on
plot([ix in]*ts,[mx mn],'g*')
plot([ix in]*ts,[0 0],'LineWidth',2)
saveas(gcf,strcat(dirname,'figure5-con.jpeg')) %exportgraphics
% 
% 
figure(6);
hold on
plot(sim_output.time_s,sim_output.f_activated,'k-','LineWidth',2)
plot(sim_output.time_s,sim_output.f_bound,'k--','LineWidth',2)
saveas(gcf,strcat(dirname,'figure6-con.jpeg')) %exportgraphics


% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output
[maf inf]= max(sim_output.muscle_force);
peak(1) = maf;
AUC(1) = sum(sim_output.muscle_force*ts);
SRX(1) = min(sim_output.M1);

% Write muscle_force to file
fid = fopen(strcat(dirname,'muscle_force'), 'a');
fprintf(fid, '%s \n', "muscle force control");
fprintf(fid, '%f %f \n', maf, inf);
fclose(fid);


% create the arrays to store the force values for control

control_force_tbl(key_indx, param_indx) = maf;

% do some plotting

for ib=1:length(dist)
pops(ib,1) = sim_output.cb_pops(inf,1,ib);
end

% commented to test without graphs
figure(7)
plot(dist,pops(:,1),'k-','LineWidth',2);
saveas(gcf,strcat(dirname,'figure7-con.jpeg')) %exportgraphics
% 
% 
figure(3);
subplot(3,1,1);
plot(sim_output.time_s,sim_output.cb_force,'k-','LineWidth',2);
ylabel('Force (N m^{-2})');
subplot(3,1,2);
plot(sim_output.time_s,sim_output.hs_length,'k-','LineWidth',2);
ylabel('Half-sarcomere length (nm)');
subplot(3,1,3);
plot(sim_output.time_s,sim_output.M1,'k-','LineWidth',2);
ylabel('Myosin in SRX (M_{OFF})');
saveas(gcf,strcat(dirname,'figure3-con.jpeg'))%exportgraphics
% 
% 
figure(4)
subplot(1,2,1);
plot(sim_output.time_s,sim_output.M2,'k-','LineWidth',2);
ylabel('Myosin available for binding (M_{ON}');
subplot(1,2,2);
N_unbound = sim_output.f_activated-sim_output.f_bound;
plot(sim_output.time_s,N_unbound,'k-','LineWidth',2);
ylabel('Thin filament available for binding (N_{UNBOUND}');
saveas(gcf,strcat(dirname,'figure4-con.jpeg'))


%Variables_1s.txt';%
protocol_file_string = 'protocol_exp_P710R2.txt';%
model_parameters_json_file_string = 'twitch_P710R_model_exp.json';%_exp
options_file_string = 'twitch_1_options.json';
model_output_file_string = '..\..\temp\twitch_P710R_output.myo';

%change the json parameters for simulation
json_struct_con_exp = loadjson('twitch_P710R_model_exp.json');%exp
if(length(parameter_var) == 2)
json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}) = (json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}))*parameter_factor;
end
if(length(parameter_var) == 3)
json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}).(parameter_var{3}) = (json_struct_con_exp.MyoSim_model.(parameter_var{1}).(parameter_var{2}).(parameter_var{3}))*parameter_factor;
end
savejson('',json_struct_con_exp,'FileName','twitch_P710R_model_exp_changed2.json')
model_parameters_json_file_string = 'twitch_P710R_model_exp_changed2.json';%exp

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

[mx, ix]=max(vel);
[mn, in]=min(vel);
con_time(2) = (ix-in)*ts
[maf inf]= max(sim_output.muscle_force);

% Write muscle_force to file
fid = fopen(strcat(dirname,'muscle_force'), 'a+');
fprintf(fid, '%s \n', "muscle force sim");
fprintf(fid, '%f %f \n', maf, inf);
fclose(fid);

% create the arrays to store the force values for control

sim_force_tbl(key_indx, param_indx) = maf;
%

peak(2) = maf
peak_ratio = peak(2)/peak(1)
AUC(2) = sum(sim_output.muscle_force*ts)
Force_Index = AUC(2)/AUC(1)
SRX(2) = min(sim_output.M1)
for ib=1:length(dist)
pops(ib,2) = sim_output.cb_pops(inf,1,ib);
end

% commented to test without graphs
figure(7)
hold on
plot(dist,pops(:,2),'r-','LineWidth',2);
ylabel('bound myosin distribution');
saveas(gcf,strcat(dirname,'figure7-sim.jpeg'))%exportgraphics
clf
% 
figure(5);
hold on
plot(sim_output.time_s,vel,'r-','LineWidth',2)
plot([ix in]*ts,[mx mn],'y*')
plot([ix in]*ts,[50,50],'LineWidth',2)
saveas(gcf,strcat(dirname,'figure5-sim.jpeg'))
clf
% 
figure(6);
hold on
plot(sim_output.time_s,sim_output.f_activated,'r-','LineWidth',2)
plot(sim_output.time_s,sim_output.f_bound,'r--','LineWidth',2)
legend('Control activated','Control bound','P710R activated','P710R bound')
saveas(gcf,strcat(dirname,'figure6-sim.jpeg')) %exportgraphics
clf

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

% commented to test without graphs
figure(3);
subplot(3,1,1);
hold on
plot(sim_output.time_s,sim_output.cb_force,'r-','LineWidth',2);
subplot(3,1,2);
hold on
plot(sim_output.time_s,sim_output.hs_length,'r-','LineWidth',2);
subplot(3,1,3);
hold on
plot(sim_output.time_s,sim_output.M1,'r-','LineWidth',2);
saveas(gcf,strcat(dirname,'figure3-sim.jpeg')) %exportgraphics
clf
% 
figure(4);
subplot(1,2,1);
hold on
plot(sim_output.time_s,sim_output.M2,'r-','LineWidth',2);
subplot(1,2,2);
hold on
N_unbound = sim_output.f_activated-sim_output.f_bound;
plot(sim_output.time_s,N_unbound,'r-','LineWidth',2);
saveas(gcf,strcat(dirname,'figure4-sim.jpeg')) %exportgraphics
clf
% 
if plot_experimental
%     t_exp = load('cell_time_trace.mat');
%     Con_F = load('Control_c4_force.mat');
%     P710R_F = load('P710R_c32_force.mat');
%     figure(3)
%     subplot(3,1,1);
%     hold on
%     plot(t_exp.cell_time_trace,Con_F.Control_c4_force,'k-');%+113.5
%     plot(t_exp.cell_time_trace,P710R_F.P710R_c32_force,'r-');%+1269.1
%     subplot(3,1,3)
%     hold on
%     plot([0.2 1],[0.60 0.6],'k:')
%     plot([0.2 1],[0.40 0.4],'r:')  
%     exportgraphics(gcf,strcat(dirname,'figure3-sim-exp.jpeg'))
%     clf
% end

 
end
pline = fgetl(parfid);

% display all the force values stored for all parameters and factors
disp(keylist)
disp(control_force_tbl)
disp(sim_force_tbl)
end
fclose(parfid);

% copy to table and save as spreadsheet
% create a table with the params and the muscle force values for each factor 
dim = (size(keylist'));
T_control  = cell2table(cell(dim(1),6), 'VariableNames', {'param', matlab.lang.makeValidName('0.1'), matlab.lang.makeValidName('0.5'), matlab.lang.makeValidName('1'), matlab.lang.makeValidName('2'), matlab.lang.makeValidName('10')});
T_sim  = cell2table(cell(dim(1),6), 'VariableNames', {'param', matlab.lang.makeValidName('0.1'), matlab.lang.makeValidName('0.5'), matlab.lang.makeValidName('1'), matlab.lang.makeValidName('2'), matlab.lang.makeValidName('10')});

% the control values of force table
T_control{:,'param'} = keylist'; % transpose keylist
T_control{:,'0.1'} = num2cell(control_force_tbl(:,1));
T_control{:,'0.5'} = num2cell(control_force_tbl(:,2));
T_control{:,'1'} = num2cell(control_force_tbl(:,3));
T_control{:,'2'} = num2cell(control_force_tbl(:,4));
T_control{:,'10'} = num2cell(control_force_tbl(:,5));

% the sim values of force table
T_sim{:,'param'} = keylist';
T_sim{:,'0.1'} = num2cell(sim_force_tbl(:,1));
T_sim{:,'0.5'} = num2cell(sim_force_tbl(:,2));
T_sim{:,'1'} = num2cell(sim_force_tbl(:,3));
T_sim{:,'2'} = num2cell(sim_force_tbl(:,4));
T_sim{:,'10'} = num2cell(sim_force_tbl(:,5));

disp(T_control)
disp(T_sim)

outputfile_control = outputdir+"param_table_control.xls"; % file name for the output control
outputfile_sim = outputdir+"param_table_sim.xls"; % file name for the output sim

writetable(T_control,outputfile_control);
writetable(T_sim,outputfile_sim);

xt = ["0.1", "0.5", "1.0", "2.0", "10.0"];
 
y = table2array(control_tbl(2,2:6));
 
x = linspace(1,5,5);
 
 
plot(x,y);
xticks([1 2 3 4 5]);
xticklabels(xt);
%set(gca,'xtick',x,'xticklabel',xt)
end
