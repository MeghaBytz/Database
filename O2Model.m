clear all
close all

%Reaction set from "Applicability of self-consistent global model for
%characterization of inductively coupled Cl2 plasma" Efremov et al

stoichiometricMatrix = [-1	0	1	0	0	0	0	-1	;
-1	1	0	0	1	0	0	-1	;
-1	1	0	1	0	0	0	-1	;
-1	0	0	1	1	0	0	-1	;
-1	2	0	0	0	0	0	-1	;
-1	1	0	0	0	0	1	-1	;
0	-1	0	1	0	0	0	-1	;
0	2	-1	0	0	0	0	-1	;
1	-1	0	0	-1	0	0	0	;
1	1	-1	0	-1	0	0	0	;
0	3	-1	0	-1	0	0	0	;
0	2	0	-1	-1	0	0	0	;
-1	1	1	-1	0	0	0	0	;
-1	0	0	0	0	1	0	-1	;
1	0	0	0	0	-1	0	-1	;
0	0	1	0	0	-1	0	-1	;
0	1	0	0	1	-1	0	-1	;
0	2	0	0	0	-1	0	0	;
1	1	0	0	-1	-1	0	-1	;
1	0	0	0	-1	-1	0	0	;
0.5	-1	0	0	0	0	0	0	;
1	0	0	0	0	-1	0	0]	;

[numberOfReactions, numberOfSpecies] = size(stoichiometricMatrix);
kRates = ones(1,numberOfReactions);
speciesConcentrations = [2 1 2 1 2];
reducedMatrix = zeros(numberOfReactions,numberOfSpecies,numberOfSpecies);
%{'a2*a5*k2 - 2*a1^2*k1'}
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
% *speciesConcentrations(2)^reducedMatrix(rxn,2,r)*speciesConcentrations(3)^reducedMatrix(rxn,3,r)...
%                 *speciesConcentrations(4)^reducedMatrix(rxn,4,r)*speciesConcentrations(5)^reducedMatrix(rxn,5,r)*speciesConcentrations(6)^reducedMatrix(rxn,6,r));
%             

%kinetic balance equations
speciesConcentrations = sym('n', [1 numberOfSpecies]);
kRates = sym('k', [1 numberOfReactions]);
eqntest = cell(numberOfSpecies,1);

for r = 1:numberOfSpecies-1
    for rxn = 1:numberOfReactions
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(speciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,['*' m1]);
            end
        end
            if ~isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(eqntest{r},[' + ' term]);
            elseif isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = term;
            end
    end
end


%mass balance eqn 

% 
% %quasineutrality condition 
% %[Cl2 Cl Cl2+ Cl+ Cl- ne]
% charges = [0 0 1 1  -1 -1];
% quasineutralityCondition = [];
% for i=1:numberOfSpecies
%     q = char(speciesConcentrations(i)*charges(i))
%     if ~isempty(quasineutralityCondition)&&~strcmp(q,'0')
%         quasineutralityCondition = strcat(quasineutralityCondition,[' + ' q]);
%     elseif isempty(quasineutralityCondition)&&~strcmp(q,'0')
%         quasineutralityCondition = q;
%     end
% end

%bohm velocity eqn




%power balance equations
