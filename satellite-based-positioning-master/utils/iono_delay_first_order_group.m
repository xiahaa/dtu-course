function d_iono = iono_delay_first_order_group(varargin)
%% iono_delay_first_order_group estimate the ionosphere delay using 
%% the first order group delay model.
% inputs:
%   el: elevation angle in degree.
%   TECU: vertical total electron content in TEC-units.
    format long;

    el = varargin{1};
    if nargin == 1
        %%from http://sesolstorm.kartverket.no/moreplots.xhtml, date 16/9/18 06:00
        TECU = 8;
    else
        TECU = varargin{2};
    end
    
    tecm=TECU*1e16;
    f1 = 1.57542e9;%% L1 frequency
    RE = 6371e3;
    hI = 350e3;
    d_ionov = 40.3/f1^2*tecm;%%omit 1e16 since tecm is normalized
    
    %% obliquity factor
    zenith = pi/2 - deg2rad(el);
    OF = (1-((RE*sin(zenith))/(RE+hI))^2)^(-1/2);
    d_iono = d_ionov * OF;
end