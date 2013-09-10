clear; clc; close all;

%%
mex compute_similarity_matrix_a.cpp;
mex modify_guiding_matrix.cpp

%%
% mex -g compute_similarity_matrix_a.cpp;
% mex -g modify_guiding_matrix.cpp