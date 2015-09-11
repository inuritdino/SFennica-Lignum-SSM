function z = lpfg_lgm_fitness_perttunen1998(X,data,model)
% Simple fitness for LIGNUM to visually correspond to the tree in Perttunen
% 1998.

[~,A] = model(X);
if(isempty(A))
    disp('Simulation had errors. Set max distance.');
    z = 10.0;
else
    % DiamBase
    tf = strcmp('DiamBase',A.colheaders);
    DiamBase = A.data(tf);
    % Height
    tf = strcmp('Height',A.colheaders);
    Height = A.data(tf);
    
    % Calculate score
    z2 = abs( (DiamBase - data(1))/data(1) );
    z3 = abs( (Height - data(2))/data(2) );
    
    z = (z2 + z3) / 2;
    
    % Plot
    figure(10);
    stem([z2 z3]);
    title(['Dist = ' num2str(z)]);
end
end