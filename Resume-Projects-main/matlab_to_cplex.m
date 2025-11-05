function [] = matlab_to_cplex( achFileName, bMin, afObjectiveCoefs, ...
                               aafAConstraintMat, afbRHS, ...
                               aafAConstraintMatEq, afbRHSEq, ...
                               afLB, afUB, achVarType )
% Purpose: Save a .lp file representing the LP problem to use with cplex
% Input Variables:
% achFileName         File name with extension .lp
% bMin                1 for a minimzation problem, 0 for a maximization problem
% afObjectiveCoefs    Column vector of objective function coefficient for each DV
% aafAConstraintMat   Matrix of DV coefficients for <= constraints
% afbRHS              Column vector of RHS values for <= constraints
% aafAConstraintMatEq Matrix of DV coefficients for equality constraints
% afbRHSEq            Column vector of RHS values for equality constraints
% afLB                Column vector of lower bound for each DV
% afUB                Column vector of upper bound for each DV
% achVarType          Column vector giving type of each DV 
%                     ('I' for integer, 'C' for real )

% In CPLEX do:
% CPLEX> read fullpath/filename.lp 
% CPLEX> opt
% CPLEX> write fullpath/filename.sol

    % Set a big number constant for inf
    fInfinity = 1e16;

    % Open the file and write preamble
    fID = fopen( achFileName, 'wt' );
    fprintf( fID, '\\Problem name: %s\n\n', achFileName );

    % Determine number of constraints and variables
    [ cConstraint, cDV ] = size( aafAConstraintMat );
    [ cConstraintEq, cDV ] = size( aafAConstraintMatEq );

    % Objective function
    if bMin
        fprintf( fID, 'Minimize\n obj: ' );
    else
        fprintf( fID, 'Maximize\n obj: ' );
    end
    fprintf( fID, '%f x%d ', afObjectiveCoefs( 1 ), 1 );
    for iDV = 2 : cDV
        if afObjectiveCoefs( iDV ) > 0
            fprintf( fID, '+ %f x%d ', afObjectiveCoefs( iDV ), iDV );
        elseif afObjectiveCoefs( iDV ) < 0
            fprintf( fID, ' %f x%d ', afObjectiveCoefs( iDV ), iDV );
        end
    end

    % Write constraints preamble
    fprintf( fID, '\nSubject To\n' );

    % Inequality <= constraints
    cConstraintWritten = 0;
    for iConstraint = 1 : cConstraint
        cConstraintWritten = cConstraintWritten + 1;
        fprintf( fID, ' c%d: %f x%d ', cConstraintWritten, aafAConstraintMat( iConstraint, 1 ), 1 );
        for iDV = 2 : cDV
            if aafAConstraintMat( iConstraint, iDV ) > 0
                fprintf( fID, '+ %f x%d ', aafAConstraintMat( iConstraint, iDV ), iDV );
            elseif aafAConstraintMat( iConstraint, iDV ) < 0
                fprintf( fID, ' %f x%d ', aafAConstraintMat( iConstraint, iDV ), iDV );
            end
        end
        if afbRHS( iConstraint ) == inf
            fprintf( fID, '<= %f\n', fInfinity );
        else
            fprintf( fID, '<= %f\n', afbRHS( iConstraint ) );
        end
    end

    % Equality = constraints
    for iConstraint = 1 : cConstraintEq
        cConstraintWritten = cConstraintWritten + 1;
        fprintf( fID, ' c%d: %f x%d ', cConstraintWritten, aafAConstraintMatEq( iConstraint, 1 ), 1 );
        for iDV = 2 : cDV
            if aafAConstraintMatEq( iConstraint, iDV ) > 0
                fprintf( fID, '+ %f x%d ', aafAConstraintMatEq( iConstraint, iDV ), iDV );
            elseif aafAConstraintMatEq( iConstraint, iDV ) < 0
                fprintf( fID, ' %f x%d ', aafAConstraintMatEq( iConstraint, iDV ), iDV );
            end
        end
        if afbRHSEq( iConstraint ) == inf
            fprintf( fID, '= %f\n', fInfinity );
        else
            fprintf( fID, '= %f\n', afbRHSEq( iConstraint ) );
        end
    end

    % Write bound preamble
    fprintf( fID, '\nBounds\n' );

    % Lower and upper bounds
    for iDV = 1 : cDV
        fprintf( fID, ' %f <= x%d <= %f\n', afLB( iDV ), iDV, afUB( iDV ) );
    end

    % Variable types
    % Binary first
    bWrittenHeader = 0;
    for iDV = 1 : cDV
        if achVarType( iDV ) == 'I' && afUB( iDV ) == 1
            if bWrittenHeader == 0
                fprintf( fID, 'Binary\n' );
                bWrittenHeader = 1;
            end
            fprintf( fID, ' x%d\n', iDV );
        end
    end
    % General integer second
    bWrittenHeader = 0;
    for iDV = 1 : cDV
        if achVarType( iDV ) == 'I' && afUB( iDV ) ~= 1
            if bWrittenHeader == 0
                fprintf( fID, 'General\n' );
                bWrittenHeader = 1;
            end
            fprintf( fID, ' x%d\n', iDV );
        end
    end

    % Write closing and close the file
    fprintf( fID, 'End\n' );
    fclose( fID );

end
