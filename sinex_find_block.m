function [lin_no, lin_pos] = sinex_find_block(fid, block_name)
% SINEX_FIND_BLOCK finds the specified block name in a SINEX file.
%   [lin_no, lin_pos] = sinex_find_block(fid, block_name) looks for the 
%   'block_name' block in the previously opened SINEX file with fid file 
%   identifier, and returns the line number and the pointer position of the 
%   input block in bytes. It returns zero for both outputs, if the
%   block_name can not be found.

% Developed by: Mohammad Ali Goudarzi (ma.goudarzi [at] gmail.com)
% Version: 1.0
% Last update: July 25, 2014

lin_no = 0;
lin_pos = 0;
block_length = length(block_name);

frewind(fid);
while ~feof(fid)
    tline = fgetl(fid);
    lin_no = lin_no + 1;
    if isequal(length(tline), block_length) && strcmpi(tline(1:block_length), block_name)
        lin_pos = ftell(fid);
        return;
    end
end

end