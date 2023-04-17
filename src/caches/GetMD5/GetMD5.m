function GetMD5(varargin)
% Calculates a 128 bit MD5 checksum for arrays or files.
%
% Usage
% -----
% :code:`digest = GetMD5(Data, Mode, Format)`
%
% Arguments
% ---------
% Data
%   input data on which the checksum is computed.
%
% Mode : :code:`char`
%   optional declaration of the type of the input data.
%
%   - 'File' : `data` is a file name as a `char`.
%   - '8Bit' : If `data` is a `char` array, only the 8 bit ASCII part is used. Then the
%   digest is the same as for a ASCII text file e.g. created by `fwrite(fid, data,
%   'uchar')`.
%   - 'Binary' : The MD5 sum is obtained for the contents of `data`. This works for
%   numerical, char and logical arrays.
%   - 'Array' : Include the class and size information of `data` in the MD5 sum. This can be
%   applied for (nested) structs, objects, cells and sparse arrays also.
%
% Format : char
%   Format of the output, default value is 'hex'.
%
%   - 'hex' : (1, 32) lowercase hexadecimal char.
%   - 'HEX' : (1, 32) uppercase hexadecimal char.
%   - 'double' : (1, 16) double vector with uint8 values.
%   - 'uint8' : (1, 16) uint8 vector.
%   - 'base64' : (1, 22) char vector, encoded in base 64 (A:Z, a:z, 0:9, +, /), unpadded.
%
% Returns
% -------
% digest
%   A 128 bit number in the specified format.
%
% Notes
% -----
% For sparse arrays, function handles, java and user-defined objects :func:`GetMD5_helper`
% is called to convert into a data format that can be handled.
%
% The C-Mex-file is compiled automatically when this function is called for the first time.
% See :file:`GetMD5.c` for details or a manual installation.
%
% References
% ----------
% Author: Jan Simon, Heidelberg, (C) 2009-2019 matlab.2010(a)n(MINUS)simon.de
%
% License: This program is derived from the RSA Data Security, Inc.
%          MD5 Message Digest Algorithm, RFC 1321, R. Rivest, April 1992
%          This implementation is published under the BSD license.
%
% See also
% --------
% Other methods for checksums can be found: :code:`CalcCRC32`, :code:`DataHash`, etc.
%
% For more checksum methods see
% `here <http://www.mathworks.com/matlabcentral/fileexchange/31272-datahash>`_.

% Dummy code, which calls the auto-compilation only: ---------------------------
% This M-function is not called, if the compiled MEX function is in the path.

persistent FirstRun
if isempty(FirstRun)
   fprintf('### %s: M-verion is running. Compiling MEX file...\n', mfilename);
   
   ok = InstallMex('GetMD5.c', 'uTest_GetMD5');
   if ok
      FirstRun = false;
   end

   try
      dummy = GetMD5_helper([]);  %#ok<NASGU>
   catch ME
      error(['JSimon:', mfilename, ':NoHelper'], ...
         '### %s: GetMD5_Helper.m: %s', mfilename, ME.message);
   end
else
   error(['JSimon:', mfilename, ':MissMEX'], ...
      'Cannot find Mex file: %s', [mfilename, '.', mexext]);
end

end
