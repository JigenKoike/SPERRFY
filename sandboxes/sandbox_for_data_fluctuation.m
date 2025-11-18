% sandbox_for_data_fluctuation.m

%% data preparation
%% Startup
clear; clc;
startup_SPERRFY;  % path setting

%% Load data
main_import_processed_data;  % import data (.mat)
main_data_preparation;
main_parameter_setting;


%% fluctuation_ConnectionBinariThreshold
connectionPValueMatrix = load()