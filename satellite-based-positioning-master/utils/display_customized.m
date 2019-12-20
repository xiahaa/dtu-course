function h1 = display_customized(varargin)
    x = varargin{1};
    y = varargin{2};
    titles = varargin{3};
    type = varargin{4};
    if type == 1
        figure;
        h1 = plot(x,y,'-o');
        grid on;
        title(titles, 'Interpreter', 'latex');
    else
        figure;
        h1 = semilogy(x,y,'-o');
        grid on;
        title(titles, 'Interpreter', 'latex');
    end
    
    if nargin > 4
        xlabels = varargin{5};
        xlabel(xlabels, 'Interpreter', 'latex');
    end
    
    if nargin > 5
        ylabels = varargin{6};
        ylabel(ylabels, 'Interpreter', 'latex');
    end
end