function [ score ] = circ_in_square( num, dim, flags, inputs )
%CIRC_IN_SQUARE objective function for the circles-in-squares optimisation
%    problem.
%    Parameters:
%        num - the number of circles in the square
%        dim - the dimensionality of the problem
%        flags - 2x1 vector of flag variables
%            flags(1) = 0 - minimisation problem, 1 - maximisation problem
%            flags(2) = 0 - handle constraint violations by repairing
%            values.  1 - handle constraint violations by returning a poor
%            score.  2 - Repair and penalise by magnitude of repair.
%        inputs - (num x dim) x 1 vector of the coordinates of the centres
%        of the circles.  coords[1:num] = x values.  coords[num+1:2*num] =
%        y values ...
%    Example SIMTOOL problem struct: (5 circles, 2D)
%       prob.id = 5
%       prob.func = @(x)circ_in_square(5, 2, [0;2],x);
%       prob.dims = 5*2;
%       prob.initgen = @(x)rand(5*2,1);
%       prob.minvals = zeros(1, 5*2);
%       prob.maxvals = ones(1, 5*2);
%       prob.optval = sqrt(2)/2;

% validate input
if size(flags, 1) ~= 2
    error('CIRC_IN_SQUARE: invalid flags');
end

if num <= 0
    error('CIRC_IN_SQUARE: number of circles must be greater than 0');
end

if dim <= 0
    error('CIRC_IN_SQUARE: number of dimensions must be greater than 0');
end

if size(inputs) ~= num * dim
    error('CIRC_IN_SQUARE: invalid coordinates matrix');
end

% Can't do much for only a single point
if num == 1
    score = 0;
    return
end

% Check for constraint violations, i.e., points with coordinate values
% outside of [0,1]
if flags(2) == 0
    % repair mode
    for i = 1:num * dim
        if inputs(i) < 0 
            inputs(i) = 0;
        elseif inputs(i) > 1
            inputs(i) = 1;
        end
    end
elseif flags(2) == 1
    % return poor score mode
    for i = 1:num * dim
        if inputs(i) < 0 || inputs(i) > 1
            % make poor score proportional to distance from valid solution
            score = sum(inputs(inputs < 0)) - sum(inputs(inputs > 1));
            % If this is a minimisation problem, flip sign of score
            if flags(1) == 0
                score = -score;
            end
            return
        end
    end
elseif flags(2) == 2
    if isempty(find(inputs < 0, 1)) && isempty(find(inputs > 1, 1))
        % No constraint violation
        penalty = 0;
    else
        % Set smallest coordinate value to at least zero
        penalty = -min(min(inputs),0);
        inputs = inputs + penalty;
        
        % Scale points so maximum coordinate is <= 1
        normfact = max(inputs);
        if(normfact > 1)
            inputs = inputs / normfact;
            % Penalise normalisation factors the further they are from 1
            penalty = penalty + normfact^2 - 1;
        end
    end
end

% Find the minimum distance between any two points
coords = reshape(inputs, num, dim);
score = min(pdist(coords));

% Subtract penalty from score
if flags(2) == 2
    score = score - penalty;
end

% If this is a minimisation problem, flip sign of score
if flags(1) == 0
    score = -score;
end
