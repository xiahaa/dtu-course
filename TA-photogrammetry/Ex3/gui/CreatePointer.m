% Taken from pixval.m:
function [pointerShape, pointerHotSpot] = CreatePointer
    pointerHotSpot = [8 8];
    pointerShape = [ ...
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
        1   1   1   1   1    1   2  NaN  2   1   1   1   1   1   1   1
        2   2   2   2   2    2  NaN NaN NaN  2   2   2   2   2   2   2
       NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
        2   2   2   2   2    2  NaN NaN NaN  2   2   2   2   2   2   2
        1   1   1   1   1    1   2  NaN  2   1   1   1   1   1   1   1
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
       NaN NaN NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN NaN NaN NaN NaN];
end