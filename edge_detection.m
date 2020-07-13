%% 
clc; clear; close all
iFish = 1;
 A = imread(['bf p' num2str(iFish),'.png']);
 roipoly(A)