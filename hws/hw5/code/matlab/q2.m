clear; clc;

N = 5;
sigma = 0.01;

[Oij, RijTrue] = GenerateObservations(N, sigma);
[Rhat, RijHat] = RotationSynchronization(Oij);
