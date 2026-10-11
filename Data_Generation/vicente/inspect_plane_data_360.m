% inspect_plane_data_360.m
%
% Lightweight schema check (whos('-file', ...) doesn't load data) on the
% plane_data_360 files, to confirm they carry the fields the
% identification pipeline (load_pole_data.m / build_glrt_model.m /
% classify_glrt_data.m) expects, and to sanity-check curve resolution
% (N) against the current curve (N=3520, per regen_poles.m header notes)
% rather than the old stale N=552 curve.

clear
clc

addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));

data_dir = "/scratch/msantana/plane_data_360";
files = dir(fullfile(data_dir, "*.mat"));

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    fprintf("\n=== %s ===\n", files(i).name);
    S = whos('-file', fpath);
    for j = 1:numel(S)
        fprintf("  %-20s size=%-20s class=%s\n", S(j).name, mat2str(S(j).size), S(j).class);
    end
end
