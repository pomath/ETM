function [solutionestimate, stnpos_idx, stnvel_idx] = sinex_read_solutionestimate(fid)
% SINEX_READ_SOLUTIONESTIMATE reads the 'solution/estimate' block of a SINEX file.
%   solutionestimate = sinex_read_solutionestimate(fid) reads the 
%   'solution/estimate' block in the previously opened SINEX file with fid 
%   file identifier, and returns it as a cell matrix 'solutionestimate'.
%   The output cell matrix has the following structure (as defined in SINEX
%   format description available at http://sopac.ucsd.edu/input/processing/
%   gamit/tables/sinex.txt):
%      1. Estimated Parameters Index
%      2. Parameter Type
%      3. Site Code
%      4. Point Code
%      5. Solution ID
%      6. Time: year
%      7. Time: day of year
%      8. Time: seconds of day
%      9. Parameter Units
%     10. Constraint Code
%     11. Parameter Estimate
%     12. Parameter Standard Deviation
% 
%  The function returns an empty cell in the case of any error.

% Developed by: Mohammad Ali Goudarzi (ma.goudarzi [at] gmail.com)
% Version: 1.2
% Last update: June 28, 2016

% Revision history:
%     1.2: (1) the textscan function was modified to recognize the empty space
%     in the time string. As a result, the time parameter is now provided
%     in three columns of 6 to 8. (2) documenting the code.

% Initializing the output parameters
solutionestimate = {};
stn_idx          = {};

% Find the line numbers for +-SOLUTION/ESTIMATE block
[start_lin, start_pos] = sinex_find_block(fid, '+SOLUTION/ESTIMATE');
end_lin                = sinex_find_block(fid, '-SOLUTION/ESTIMATE');
if le(start_lin, 0) || le(end_lin, 0) || le(end_lin, start_lin)
    return;
end

% Put the curser at the begining of the +SOLUTION/ESTIMATE block
fseek_status = fseek(fid, start_pos, 'bof');
if ~isequal(fseek_status, 0)
    return;
end

tline = fgetl(fid);
if tline == -1
    return;
end

% Read the +SOLUTION/ESTIMATE block
solutionestimate = textscan(fid, '%u %s %s %s %s %f %f %f %s %u %f %f', ...
    end_lin - start_lin - 2, ...
    'delimiter', {' ', ':'}, 'MultipleDelimsAsOne', 1);

% Make the index of the stations in the +SOLUTION/ESTIMATE block. Other
% parameters than stations' coordinates might be there.
indx = strcmpi('STAX', solutionestimate{2});
indy = strcmpi('STAY', solutionestimate{2});
indz = strcmpi('STAZ', solutionestimate{2});
stnpos_idx = logical(sum([indx, indy, indz], 2));

indx = strcmpi('VELX', solutionestimate{2});
indy = strcmpi('VELY', solutionestimate{2});
indz = strcmpi('VELZ', solutionestimate{2});
stnvel_idx = logical(sum([indx, indy, indz], 2));
end