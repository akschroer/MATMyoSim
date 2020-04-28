function demo_ramp_2
% Function illustrates how to run a simulation of a single-half-sarcomere
% connected in series with a spring, activated at a constant pCa
% Myosin heads cycle but do not generate active force

% Variables
protocol_file_string = '..\..\temp\ramp_2_protocol.txt';
model_parameters_json_file_string = 'ramp_2_parameters.json';
options_file_string = 'ramp_2_options.json';
model_output_file_string = '..\..\temp\ramp_2_output.myo';

% Make sure the path allows us to find the right files
addpath('..\..\..\generate_protocols','..\..\..\..\code');

% Generate a new protocol
generate_triangle_protocol( ...
    'output_file_string', protocol_file_string);

% Run a simulation
sim_output = simulation_driver( ...
    'simulation_protocol_file_string', protocol_file_string, ...
    'model_json_file_string', model_parameters_json_file_string, ...
    'options_file_string', options_file_string, ...
    'output_file_string', model_output_file_string);

% Load it back up and display to show how that can be done
sim = load(model_output_file_string,'-mat')
sim_output = sim.sim_output

figure(2);
clf;
subplot(2,1,1);
plot(sim_output.time_s,sim_output.muscle_force,'b-');
ylabel('Force (N m^{-2})');
subplot(2,1,2);
plot(sim_output.time_s,sim_output.hs_length,'b-');
ylabel('Half-sarcomere length (nm)');