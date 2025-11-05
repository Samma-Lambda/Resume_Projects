function [ afDV, fObjValue ] = cplex_to_matlab( achFileName )
% Purpose: Load cplex results from a .sol file
% Input Variables:
% achFileName CPLEX xml output file name (generally with extension .sol)
% Output Variables:
% afDV        Vector of optimal value for each DV
% fObjValue   Optimal value of objective function

% In CPLEX do:
% CPLEX> read fullpath/filename.lp 
% CPLEX> opt
% CPLEX> write fullpath/filename.sol

    % Open the file
    fID = fopen( achFileName );
    
    % Read the file one line at a time into a string variable for the line
    achLine = fgetl( fID );
    while achLine ~= -1
        % If the line contains the objective function value
        nIndex = strfind( achLine, 'objectiveValue="' );
        if nIndex > 0
            achFrontCropped = achLine( nIndex + 16 : length( achLine ) );
            anIndices = strfind( achFrontCropped, '"' );
            achObjValue = achFrontCropped( 1 : min( anIndices ) - 1 );
            fObjValue = strread( achObjValue, '%f' );
        end
        % If the line contains variable name
        nIndex = strfind( achLine, '<variable name="x' );
        if nIndex > 0
            % Determine the DV number as an integer
            achFrontCropped = achLine( nIndex + 17 : length( achLine ) );
            anIndices = strfind( achFrontCropped, '"' );
            achDVNum = achFrontCropped( 1 : min( anIndices ) - 1 );
            nDVNum = strread( achDVNum, '%d' );
            % Determine the DV value
            nIndex = strfind( achFrontCropped, 'value="' );
            achMoreCropped = achFrontCropped( nIndex + 7 : length( achFrontCropped ) ) ;
            anIndices = strfind( achMoreCropped, '"' );
            achDVValue = achMoreCropped( 1 : min( anIndices ) - 1 );
            fDVValue = strread( achDVValue, '%f' );
            % Fill in value for current DV in output vector
            afDV( nDVNum ) = fDVValue;
        end
        % Read the next line
        achLine = fgetl( fID );
    end
 
    % Close the file
    fclose( fID );

end
