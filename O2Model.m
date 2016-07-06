clear all
close all

%Reaction set from "Applicability of self-consistent global model for
%characterization of inductively coupled Cl2 plasma" Efremov et al
charges = [0 0 1 1  -1 0 0 -1];
gasInput = [1 0 0 0 0 0 0 0 0]
stoichiometricMatrix = [-1	0	1	0	0	0	0	-1	0	;
-1	1	0	0	1	0	0	-1	0	;
-1	1	0	1	0	0	0	-1	0	;
-1	0	0	1	1	0	0	-1	0	;
-1	2	0	0	0	0	0	-1	0	;
-1	1	0	0	0	0	1	-1	0	;
0	-1	0	1	0	0	0	-1	0	;
0	2	-1	0	0	0	0	-1	0	;
1	-1	0	0	-1	0	0	0	0	;
-1	1	1	-1	0	0	0	0	0	;
-1	0	0	0	0	1	0	-1	0	;
1	0	0	0	0	-1	0	-1	0	;
0	0	1	0	0	-1	0	-1	0	;
0	1	0	0	1	-1	0	-1	0	;
0	2	0	0	0	-1	0	0	0	;
1	1	0	0	-1	-1	0	0	0	;
0.5	-1	0	0	0	0	0	0	-1	;
1	0	0	0	0	-1	0	0	-1	;
1	0	-1	0	0	0	0	0	-1	;
0	1	0	-1	0	0	0	0	-1	;
1	1	-1	0	-1	0	0	0	0	;
0	3	-1	0	-1	0	0	0	0	;
0	2	0	-1	-1	0	0	0	0	];



numberOfChargedSpecies = nnz(charges)-1; %subtract one because we are not counting electrons
numberOfNeutrals = length(charges)-nnz(charges);
numberOfPositve = length(find(charges>0));
numberOfNegative = length(find(charges<0))-1;%not counting electrons
numberOfInputs = nnz(gasInput);


[numberOfReactions, numberOfSpecies] = size(stoichiometricMatrix);
kRates = ones(1,numberOfReactions);
speciesConcentrations = [2 1 2 1 2];
reducedMatrix = zeros(numberOfReactions,numberOfSpecies,numberOfSpecies);

%Flag recombination rxns
recReactionIndices = [];
wallReactionIndices = [];
for k = 1:numberOfReactions
    nonZeroElements = find(stoichiometricMatrix(k,1:end-2)<0); %not including ne or wall
    f = 0;
    flagged = 0;
    while flagged == 0 && f < length(nonZeroElements)
        f = f + 1;
        if charges(nonZeroElements(f))>0
            for g = 1:length(nonZeroElements)
                if charges(nonZeroElements(g))<0; 
                    recReactionIndices = [recReactionIndices k];
                    flagged = 1;
                end
            end
        end
    end
end

for k = 1:numberOfReactions
        if stoichiometricMatrix(k,end)==-1
            wallReactionIndices = [wallReactionIndices k];
        end
end

numberOfWallReactions = length(wallReactionIndices);
numberOfRecombinationReactions = length(recReactionIndices);
for i = 1:length(recReactionIndices)
    recMatrix(i,:) = stoichiometricMatrix(recReactionIndices(i),:);
end
for i = 1:length(wallReactionIndices)
    wallMatrix(i,:) = stoichiometricMatrix(wallReactionIndices(i),:);
end

stoichiometricMatrix(recReactionIndices,:) = [];
stoichiometricMatrix(wallReactionIndices,:) = [];
%add matrix back in
stoichiometricMatrix = [stoichiometricMatrix;wallMatrix];
stoichiometricMatrix = [stoichiometricMatrix;recMatrix];
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



%kinetic balance equations
volumeSpeciesConcentrations = sym('n', [1 numberOfSpecies]);
peakSpeciesConcentrations = sym('np', [1 numberOfSpecies]);
gasInputs = sym('Q', [1 numberOfInputs]);
kRates = sym('k', [1 numberOfReactions]);
eqntest = cell(numberOfSpecies,1);

for r = 1:numberOfSpecies-2 %not doing reactions for walls or ne
    for rxn = 1:numberOfReactions-(numberOfRecombinationReactions + numberOfWallReactions)
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(volumeSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(eqntest{r},[' + ' term]);
            elseif isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = term;
            end
    end
    if charges(r) == 0
        eqntest{r} = strcat(eqntest{r},[' - Kpump*' char(volumeSpeciesConcentrations(r))]);
    end
end

%add in Krec terms
for r = 1:numberOfSpecies-1
    for rxn = numberOfReactions-numberOfRecombinationReactions:numberOfReactions
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(peakSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(eqntest{r},[' + ' term ' * Vrec/V']);
            elseif isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(term,'*Vrec/V');
            end
    end

end

%add in wall reactions
for r = 1:numberOfSpecies-1
    for rxn = numberOfReactions-(numberOfRecombinationReactions + numberOfWallReactions):numberOfReactions
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(peakSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(eqntest{r},[' + ' term]);
            elseif isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = term;
            end
    end

end

% add molecular flow rates
for r = 1:numberOfSpecies-1
    if gasInput(r) ==1
      term = char(gasInputs(r));
            if ~isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(eqntest{r},[' + ' term '/V']);
            elseif isempty(eqntest{r})&&~strcmp(term,'0')
                eqntest{r} = strcat(term, '/V');
            end
    end

end


%Power balance eqn
%Ec_O2*Kiz1*nO2*ne0 - Ec_O*Kiz2*nObar*ne0 - Eei_O*Kion*nOplus - Eei_O2*Kion*nO2plus; 

collReactionIndices = [];
for k = 1:numberOfReactions-2
    for i = 1:numberOfNeutrals
            if stoichiometricMatrix(k,i)<0
                if stoichiometricMatrix(k,i+numberOfNeutrals)>0 && nnz(stoichiometricMatrix(k,numberOfNeutral+1:numberOfNeutrals+numberOfPositive))==1
                    collReactionIndices = [collReactionIndices k];
                end
            end
    end 
end





























myValues = ['k2 '; 'k3 '; 'k4 '; 'k5 '; 'k6 '; 'k7 '; 'k8 '; 'k9 '; 'k10'; 'k11'; 'k12'; 'k13'; 'k14'; 'k15'; 'k16'; 'k17'; 'k18'; 'k19'; 'k20'; 'k21'; 'k22'; 'k23';'k1 ']
cellMyValues = cellstr(myValues);
paperValues = ['Katt '; 'Kiz4 '; 'Kiz3 '; 'Kdiss'; 'Kdiss'; 'Kiz2 '; 'Kei  '; 'Kdet '; 'Kch  '; 'Kex  '; 'Kdeex'; 'Kizm '; 'Kattm'; 'Kdism'; 'Krec4'; 'vO   '; 'vO2* '; 'vO2+ '; 'vO+  '; 'Krec '; 'Krec2'; 'Krec3';'Kiz1 '];
cellPaperValues = cellstr(paperValues);
mySpecies = ['n1 '; 'n2 '; 'n3 ';'n4 '; 'n5 '; 'n6 '; 'n7 '; 'n8 ';'np3';'np4';'np5'];
cellMySpecies = cellstr(mySpecies);
paperSpecies = ['nO2  '; 'nO   '; 'nO2+ ';'nO+  '; 'nO-  '; 'nO2* '; 'nO*  '; 'ne   ';'nO2+p';'nO+p ';'nO-p '];
cellPaperSpecies = cellstr(paperSpecies);
% for j = 1:numberOfSpecies-1
%     for i=1:numberOfReactions
%         eqntest{j} =strrep(eqntest{j},cellMyValues(i),cellPaperValues(i));
%     end
% end

% for j = 1:numberOfSpecies-1
%     for i=1:length(mySpecies)
%         eqntest{j} =strrep(eqntest{j},cellMySpecies(i),cellPaperSpecies(i));
%     end
% end


%mass balance eqn 

% 
% %quasineutrality condition 
% %[Cl2 Cl Cl2+ Cl+ Cl- ne]
% 
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
