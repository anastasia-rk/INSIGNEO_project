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
%% Results folder
FigFolder   = 'Results/Figures/';
TikzFolder  = 'Results/Tikzes/';
TabFolder   = 'Results/Tables/';
%% Figure and interpreter setup
my_init;
saveFormat = '-dpng'; % '-depsc' %
visFlag = 'On'; %% or 'Off'
disFlag = true; % set to false to prevent from displaying updates in work space
