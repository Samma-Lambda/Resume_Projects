function [ afOptDV, fMinObj ] = cplex_solve( ...
    bMin, afObjectiveCoefs, aafAConstraintMat, afbRHS, ...
    aafAConstraintMatEq, afbRHSEq, afLB, afUB, VarType, ...
    achFileNameBase, achDirectory, bShowOutput, bWindows )
% Solve a linear programming problem by running CPLEX via a system call
%
% Input Variables
% bMin                 1 for a minimization problem, 0 for maximization
% afObjectiveCoefs     n x 1 vector of objective function coefficients
% aafAConstraintMat    m x n matrix of DV coefficients for m <= constraints
% afbRHS               m x 1 RHS vector for m <= constraints
% aafAConstraintMatEq  p x n matrix of DV coefficients for p = constraints
% afbRHSEq             p x 1 RHS vector for p = constraints
% afLB                 n x 1 vector of the lower bound for each DV
% afUB                 n x 1 vector of the upper bound for each DV
% VarType              n x 1 vector of DV type ('I' = integer, 'C' = real)
% achFileNameBase      Name for output files (without an extension)
% achDirectory         Directory for output files (with full path)
% bShowOutput          1 to show CPLEX output, 0 to not
% bWindows             1 if using Windows, 0 if using Linux/Mac
%
% Output Variables
% afOptDV              n x 1 vector of optimal DV values
% fMinObj              Optimal objective function value

    % Set parameters for Windows vs Linux
    if bWindows
        chSlash = '\';
        achRemove = 'del ';
    else
        chSlash = '/';
        achRemove = 'rm -f ';
    end

    % Create CPLEX LP file
    matlab_to_cplex( ...
        sprintf( '%s%c%s.lp', achDirectory, chSlash, achFileNameBase ), ...
        bMin, afObjectiveCoefs, aafAConstraintMat, afbRHS, ...
        aafAConstraintMatEq, afbRHSEq, afLB, afUB, VarType );

    % Create CPLEX command file
    fID = fopen( sprintf( '%s%ccplex_commandfile', achDirectory, chSlash ), ...
                 'wt' );
    fprintf( fID, 'read %s%c%s.lp\n', achDirectory, chSlash, achFileNameBase );
    fprintf( fID, 'opt\n' );
    fprintf( fID, 'write %s%c%s.sol\n', achDirectory, chSlash, achFileNameBase );
    fclose( fID );

    % Remove old solution file
    system( sprintf( '%s %s%c%s.sol', achRemove, achDirectory, chSlash, ...
            achFileNameBase ) );

    % Call CPLEX to run .lp file to make .sol file
    achSystemCmd = sprintf( 'cplex -f %s%ccplex_commandfile', ...
                            achDirectory, chSlash );
    [ ~, achOutput ] = system( achSystemCmd );
    if bShowOutput
        fprintf( '---------- BEGIN CPLEX OUTPUT ----------' );
        for iIndex = 1 : length( achOutput )
            if achOutput( iIndex ) == '%'
                achOutput( iIndex ) = ' ';
            end
        end
        fprintf( achOutput );
        fprintf( '\n---------- END CPLEX OUTPUT ----------\n' );
    end

    % Read results from .sol file
    [ afOptDV, fMinObj ] = ...
        cplex_to_matlab( sprintf( '%s%c%s.sol', ...
                         achDirectory, chSlash,achFileNameBase ) );

end
