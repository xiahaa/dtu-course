function deci = deg2deci(varargin)
%   deci = deg2deci(varargin) does the conversion from degree/minute/second
%   to decimal number.
%   varargin: 
%       1: do degree conversion.
%       2: 1-degree, 2-minutes.
%       3: 1-degree, 2-minutes, 3-seconds.
%   Author: xiahaa@space.dtu.dk

    numarg = length(varargin);
    if numarg == 1
        deci = varargin{1};
    elseif numarg == 2
        deci =  varargin{1} +  varargin{2}/60;
    elseif numarg == 3
        deci =  varargin{1} +  varargin{2}/60 + varargin{3}/3600;
    else
        disp('Error Format, Please enter degree, minute, second!');
    end
end