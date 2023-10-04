function [y] = TrialFunction(x,inputNN,varargin)
    % trial function for the given IVP 
    % initial condition at x=a is y(a)= ya
    % default values are y = 0, ya = 0
    % initial condition can be given as 'a',0.5
    
    % TODO: exceptions
    
    ya = 0;
    a = 0;
    for i = 1:2:length(varargin) 
        if varargin{i} == 'a'
            a = varargin{i+1};
        elseif varargin{i} == 'ya'
            ya = varargin{i+1};
        end
    end
    y = ya + (x-a) * inputNN;
end

