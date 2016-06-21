clear all
close all



stoichiometricMatrix = [-2 0 0 0 1; 1 -1 1 1 -1];
numberOfReactions = 2;
numberOfSpecies = 5;
kRates = [1 1];
speciesConcentrations = [2 1 2 1 2];
reducedMatrix = zeros(numberOfReactions,numberOfSpecies,numberOfSpecies);

for i = 1:numberOfReactions 
    for j = 1:numberOfSpecies 
            if stoichiometricMatrix(i,j)<0
                 for k = 1:numberOfSpecies 
                    if stoichiometricMatrix(i,k)<0
                         reducedMatrix(i,k,j) = abs(stoichiometricMatrix(i,k));
                    else
                         reducedMatrix(i,k,j) = 0;
                    end
                 end
             end
            if stoichiometricMatrix(i,j)>0
                for q = 1:numberOfSpecies 
                    if stoichiometricMatrix(i,q)<0
                        reducedMatrix(i,q,j) = abs(stoichiometricMatrix(i,q));
                    else
                        reducedMatrix(i,q,j) = 0;
                    end
                end
            end
            
        
    end
end

speciesConcentrations = sym('a', [1 numberOfSpecies])
kRates = sym('k', [1 numberOfReactions])
for r = 1:numberOfSpecies
            eqn{r} = char(stoichiometricMatrix(1,r)*kRates(1)*speciesConcentrations(1)^reducedMatrix(1,1,r)*speciesConcentrations(2)^reducedMatrix(1,2,r)*speciesConcentrations(3)^reducedMatrix(1,3,r)...
                *speciesConcentrations(4)^reducedMatrix(1,4,r)*speciesConcentrations(5)^reducedMatrix(1,5,r)+ stoichiometricMatrix(2,r)*kRates(2)*speciesConcentrations(1)^reducedMatrix(2,1,r)...
                *speciesConcentrations(2)^reducedMatrix(2,2,r)*speciesConcentrations(3)^reducedMatrix(2,3,r)*speciesConcentrations(4)^reducedMatrix(2,4,r)...
                *speciesConcentrations(5)^reducedMatrix(2,5,r));
end
%create 3D matrix
%how to write out all eqns
%steady state eqn


