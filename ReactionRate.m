clear all
close all

stoichiometricMatrix = [-2 0 0 0 1; 1 -1 1 1 -1]
zerosMatrix = [1 0 0 0 1; 1 1 1 1 1]
numberOfReactions = 2;
numberOfSpecies = 5;
k = [1 1];
speciesConcentrations = [ 2 1 1 1 1; 2 1 1 1 1];
for j = 5:5%numberOfSpecies
    %for i = 1:1%numberOfReactions
        test = k*speciesConcentrations.^abs(stoichiometricMatrix).*stoichiometricMatrix;
    %end
end