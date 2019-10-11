%-------------------------------------------------------------------------%
%
% Usage: run Demo.m
%
% Goal: script that gives an example for the computation of the 
%       RVSK-PU interpolant
%
% Calls on: PU_RBF_VSK_RATIONAL_DAGC (script that performs partition of 
%           unity rational interpolation via VSKs and DACG)
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%

clear all
close all
clc
warning off
addpath('RvskPuDacg') % The routine that computes the RVSK-PU via DACG
addpath('Data_Structure') % The routines that organize points on patches
addpath('Distance_Matrix') % The routines that construct distance matrices 
% for CSRBFs and VSKs
addpath('Dacg') % The routines that compute eigenvectors via DACG
addpath('Data') % The data used in this example
% Define the nodes
load('data_demo');
M = 2; % Space dimension   
N = 33^M; % Number of data
neval = 40; % Parameter for evaluation points
f = @(x) (tanh(9*(x(:,2)-x(:,1)))+1)/(tanh(9)+1); % The test function
rbf = @(ep, r) exp(-(ep*r).^2); % The kernel function
wf = @(e,r) r.^4.*(5*spones(r)-4*r); % The PU weights
npu = ceil(sqrt(N)/(4)); % Parameter for PU centres
h = 9; % Parameter for the scale function (VSKs are used)

disp('-------------   RVSK-PU rational interpolation   -----------------')
[epoints Pf] = PU_RBF_VSK_RATIONAL_DAGC(M,dsites,neval,npu,rbf,wf,f,h);