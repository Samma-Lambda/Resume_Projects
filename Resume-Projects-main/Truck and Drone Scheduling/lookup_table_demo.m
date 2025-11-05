function [ ] = lookup_table_demo( )
% This program demonstrates how to build DV forward and reverse lookup tables
% converting between multiple subscript indexing and Matlab's single-indexed
% long list of DVs

    % Set number of actual customers, trucks, and drones
    cCustActual = 5;
    cTrucks = 3;
    cDrones = 2;

    % Init: lookup tables
    [ aaanXDVLookup, aaaanYDVLookup, aanSTDVLookup, aaaanSUDVLookup, ...
      aaaanRDVLookup, aaanPDVLookup, aanUDVLookup, aaanTUADVLookup, ...
      aaanTULDVLookup, aanTTADVLookup, aanTTLDVLookup, aanWUDVLookup, ...
      aanQDVLookup, nTDVLookup, astrDVReverseLookup, aanDVReverseLookup, ...
      cDV ] = ...
        create_dv_lookup_tables( cCustActual, cTrucks, cDrones );

end

% ------------------------------------------------------------------------------------------------- %

function [ anVL ] = get_set_VL( cCustActual )
% This function returns a vector containing the elements of the set V_L using 1-indexing
% 0-indexed version: V_L = { 0, 1, 2, ..., n } = set of potential trip starting customers
% 1-indexed version: V_L = { 1, 2, 3, ..., n + 1 } = set of potential trip starting customers
    anVL = 1 : cCustActual + 1;
end

% ------------------------------------------------------------------------------------------------- %

function [ anVR ] = get_set_VR( cCustActual )
% This function returns a vector containing the elements of the set V_R using 1-indexing
% 0-indexed version: V_R = { 1, 2, ..., n + 1 } = set of potential trip ending customers
% 1-indexed version: V_R = { 2, 3, ..., n + 2 } = set of potential trip ending customers
    anVR = 2 : cCustActual + 2;
end

% ------------------------------------------------------------------------------------------------- %

function [ anC ] = get_set_C( cCustActual )
% This function returns a vector containing the elements of the set C using 1-indexing
% 0-indexed version: C = { 1, 2, ..., n } = set of customers
% 1-indexed version: C = { 2, 3, ..., n + 1 } = set of customers
    anC = 2 : cCustActual + 1;
end

% ------------------------------------------------------------------------------------------------- %

function [ aaanXDVLookup, aaaanYDVLookup, aanSTDVLookup, aaaanSUDVLookup, ...
      aaaanRDVLookup, aaanPDVLookup, aanUDVLookup, aaanTUADVLookup, ...
      aaanTULDVLookup, aanTTADVLookup, aanTTLDVLookup, aanWUDVLookup, ...
      aanQDVLookup, nTDVLookup, ...
      astrDVReverseLookup, aanDVReverseLookup, cDV  ] = ...
        create_dv_lookup_tables( cCustActual, cTrucks, cDrones )
% This function creates a lookup table for each of the DV with multiple sub/superscripts, 
% giving their DV number in a single list of DVs with one subscript.
% It also creates the reverse lookup tables which take a single DV subscript and provide the
% DV letters name and associated subscripts.

    % Initialize DV counter and customer sets
    iDV = 0;
    anVL = get_set_VL( cCustActual );
    anVR = get_set_VR( cCustActual );
    anC = get_set_C( cCustActual );

    % Set up X DV lookup table
    % x_{ijg}      (binary)      % 1 iff truck g travels from i to j; g \in \mathcal{R}, i \in V_L, j \in V_R \setminus { i }
    aaanXDVLookup = zeros( length( anVL ), length( anVR ), cTrucks );
    for iTruck = 1 : cTrucks
        for iCustI = anVL
            for iCustJ = anVR
                if iCustI ~= iCustJ
                    iDV = iDV + 1;
                    aaanXDVLookup( iCustI, iCustJ, iTruck) = iDV;
                    fprintf( 'x_{ %2d, %2d, %2d } is DV %4d\n', iCustI, iCustJ, iTruck, iDV );
                    astrDVReverseLookup{ iDV } = "x";
                    aanDVReverseLookup( iDV, : ) = [ iCustI, iCustJ, iTruck, -1 ];
                end
            end
        end
    end

    % Set up Y DV lookup table
    % y_{ijog}     (binary)      % 1 iff drone o of truck g travels from i to j; o \in \mho, g \in \mathcal{R}, i \in V_L, j \in V_R \setminus { i }
    aaaanYDVLookup = zeros( length( anVL ), length( anVR ), cDrones, cTrucks );
    for iTruck = 1 : cTrucks
        for iDrone = 1 : cDrones
            for iCustI = anVL
                for iCustJ = anVR
                    if iCustI ~= iCustJ
                        iDV = iDV + 1;
                        aaaanYDVLookup( iCustI, iCustJ, iDrone, iTruck ) = iDV;
                        fprintf( 'y_{ %2d, %2d, %2d, %2d } is DV %4d\n', ...
                                 iCustI, iCustJ, iDrone, iTruck, iDV );
                        achDVReverseLookup{ iDV } = "y";
                        aanDVReverseLookup( iDV, : ) = [ iCustI, iCustJ, iDrone, iTruck ];
                    end
                end
            end
        end
    end

    % Set up sT DV lookup table
    % sT_{ig}      (binary)      % 1 iff cust i is served by truck g; i \in C, g \in \mathcal{R}
    aanSTDVLookup = zeros( length( anC ), cTrucks );
    for iTruck = 1 : cTrucks
        for iCust = anC
            iDV = iDV + 1;
            aanSTDVLookup( iCust, iTruck ) = iDV;
            fprintf( 'sT_{ %2d, %2d } is DV %4d\n', iTruck, iCust, iDV );
            achDVReverseLookup{ iDV } = "sT";
            aanDVReverseLookup( iDV, : ) = [ iCust, iTruck, -1, -1 ];
        end
    end

    % ... more in here for other DVs ...
    aaaanSUDVLookup = [ ];
    aaaanRDVLookup = [ ];
    aaanPDVLookup = [ ];
    aanUDVLookup = [ ];
    aaanTUADVLookup = [ ];
    aaanTULDVLookup = [ ];
    aanTTADVLookup = [ ];
    aanTTLDVLookup = [ ];
    aanWUDVLookup = [ ];
    aanQDVLookup = [ ];
    nTDVLookup = [ ];
    % ... more in here for other DVs ...

    % Set final count of total DVs
    cDV = iDV;

end

% ------------------------------------------------------------------------------------------------- %
