function [p,e,t] = get_mesh_data(str)

% [p,e,t] = get_mesh_data(str)
% Take the mesh data from PDE_tool export (p,e,t) in the right subdir.
% meshes/str (where str is defined in main 
% and could be 'narrow_channel', 'circle', etc.)

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project

str_meshes = 'meshes/';
str = fullfile(str_meshes,str);
p = xlsread(fullfile(str,'/p.xlsx'));
e = xlsread(fullfile(str,'/e.xlsx'));
t = xlsread(fullfile(str,'/t.xlsx'));

end

