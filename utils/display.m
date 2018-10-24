function h1 = display(varargin)
    x = varargin{1};
    y = varargin{2};
    titles = varargin{3};
    type = varargin{4};
    if type == 1
        h1 = figure;
        plot(x,y, 'Parent', h1, 'Interpreter', 'latex');
        grid on;
        title(titles, 'Parent', h1, 'Interpreter', 'latex');
    else
        h1 = figure;
        semilogy(x,y,'Parent', h1, 'Interpreter', 'latex');
        grid on;
        title(titles, 'Parent', h1, 'Interpreter', 'latex');
    end
    
    if nargin > 4
        xlabels = varargin{5};
        xlabel(xlabels,'Parent', h1, 'Interpreter', 'latex');
    end
    
    if nargin > 5
        ylabels = varargin{6};
        ylabel(ylabels,'Parent', h1, 'Interpreter', 'latex');
    end
end