%% 9/10/14 Single Leg Example

clc;
clear;
close all;


% speicify a robot model to verify the concept
% define a 6-DOF leg
model = SingleLegYan();
q = zeros(model.NJ,1);
q(model.J_LHP) = -pi/4;
q(model.J_LKP) = pi/2;%pi/6;
q(model.J_LAP) = -pi/4;%-pi/12;
q_FB = zeros(6,1);
q_FB(1:3) = [0 0 0]; % initial position of the floating base
q = [q_FB; q];
qd = zeros(model.N_joint,1);
qdd = zeros(model.N_joint,1);

plot_config = 1;
model = model.UpdateModel(q, plot_config);

M = model.computeH()
CandG = model.computeCandG(qd)