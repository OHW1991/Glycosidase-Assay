%% Pipeline.m
% Main script for old (20251125) + new (20251130) analyses.
% Assumes folder layout:
%   RawData/
%       20251125/R8L/<conds>/FoV*/c1..c3.tif
%       20251125/R-L/<conds>/FoV*/c1..c3.tif
%       20251130/<samples>/FoV*/c1..c3.tif

close all; clear; clc;

rawRoot    = fullfile(pwd, 'RawData');
eventsRoot = fullfile(pwd, 'EventsAnalysis');

map = mapRawData(rawRoot);

%analyzeFoVEvents(map, eventsRoot);
%analyzeFoVDistances_LAMP1('RawData');

%plotRGandRB_counts_R8L_SNAplus('EventsAnalysis');
plotCoGDistanceSummaries_LAMP1('AnalysisDistance');