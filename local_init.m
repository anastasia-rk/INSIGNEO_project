clear; clc; close all;
%% Add mylibs to and subdirectories to path in local directory
% Requires that library is located in the PRESENT folder.
LibFolder = 'Functions';                                                    % name of library
pathCell = regexp(path, pathsep, 'split')';                                 % array of directories in Matlab path
index_path = find(cell2mat(strfind(pathCell,LibFolder)));                   % find whether path contains mylibs
if isempty(index_path)
    addpath(LibFolder)
    newPath = true;
else
    myPath = pathCell(index_path);                                          % see all subridectories in the path
    newPath = false;
end
%% Go through the subdirectories and add them to path
f = dir(LibFolder);
git = 1;                                                                    % git offset
N = length(f);
for iDir=git:N
    if f(iDir).isdir
       if newPath || ~any(cell2mat(strfind(myPath,f(iDir).name)))
        addpath([LibFolder,'/',f(iDir).name]);
       end
    end
end
%% Figure and interpreter setup
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultfigurecolor',[1 1 1]);
set(0,'defaulttextinterpreter','latex');  
set(0, 'defaultAxesTickLabelInterpreter','latex');  
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultColorbarTickLabelInterpreter','latex');
set(0, 'DefaultAxesFontWeight', 'normal','DefaultAxesFontSize', 14);
saveFormat = '-dpng'; % '-depsc' %
visFlag = 'On'; %% or 'Off'
disFlag = true; % set to false to prevent from displaying updates in work space

