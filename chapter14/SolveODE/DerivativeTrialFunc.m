function [yp] = DerivativeTrialFunc(x,inputNN,inputNNp,varargin)
    % derivative of trial function for the given IVP w.r.t x
    % initial condition at x=a is y(a)= ya
    % default values are y = 0
    % initial condition can be given as 'a',0.5
    
    % TODO: exceptions
    
    a = 0;
    for i = 1:2:length(varargin) 
        if varargin{i} == 'a'
            a = varargin{i+1};
        end
    end
    yp = inputNN + (x-a)*inputNNp;
end
