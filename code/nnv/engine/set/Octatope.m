classdef Octatope
    % Octatope defined by x = c + a[1]*v[1] + a[2]*v[2] + ... + a[n]*v[n]
    %                       = V * b, V = [c v[1] v[2] ... v[n]]
    %                                b = [1 a[1] a[2] ... a[n]]^T
    %                       where C*a <= d, constraints on a[i] and the constraints form a UCS


    properties
        V = []; % basic matrix
        C = []; % constraint matrix
        d = []; % constraint vector
        dim = 0; % dimension of octatope
        nVar = 0; % number of variable in the constraints
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
        state_lb = []; % lower bound of state variables
        state_ub = []; % upper bound of state variables
    end

    % constructor and main methods
    methods

        % constructor
        function obj = Octatope(varargin)
            % @V: basic matrix
            % @C: constraint matrix
            % @d: constraint vector

            switch nargin

                case 7
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    state_lb = varargin{6};
                    state_ub = varargin{7};

                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    [n1, m1] = size(pred_lb);
                    [n2, m2] = size(pred_ub);
                    [n3, m3] = size(state_lb);
                    [n4, m4] = size(state_ub);

                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column')
                    end

                    if m1 ~= 1 || m2 ~= 1
                        error('predicate lower or upper bounds vector should have one column');
                    end

                    if n1 ~= n2 || n1 ~= mC
                        error('Inconsistency between number of predicate variables and predicatre lower or upper bounds vector');
                    end

                    if n3 ~= nV || n4 ~= nV
                        error('Inconsistenct dimension between lower bound and upper bound vector of state variables and matrix V');
                    end

                    if m3 ~= 1 || m4 ~= 1
                        error('Invalid lower bound or upper bound vector of state variables');
                    end

                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                    obj.state_lb = state_lb;
                    obj.state_ub = state_ub;
            end
        end

        % affine mapping of an octatope S = Wx + b;
        function O = affineMap(obj, W, b)
            % @W: mapping matrix
            % @b: mapping vector
            % @O: new octatope
            
            if size(W, 2) ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the octatope');
            end

            if ~isempty(b)
                if size(b, 1) ~= size(W, 1)
                    error('Inconsistency between the mapping vec and mapping matrix');
                end
                if size(b, 2) ~= 1
                    error('Mapping vector should have one column');
                end
                newV = W * obj.V;
                newV(:, 1) = newV(:, 1) + b;

            else
                newV = W * obj.V;
            end

            O = Octatope(newV, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, obj.state_lb, obj.state_ub);
        end

        % intersection with a half space: H(x) := Hx <= g
        function O = intersectHalfSpace(obj, H, g)
            % @H: HalfSpace matrix
            % @g: HalfSpace vector
            % @O: new octatope with more constraints

            [nH, mH] = size(H);
            [ng, mg] = size(g);

            if mg ~= 1
                error('Halfspace vector should have one column');
            end

            if nH ~= ng
                error('Inconsistent dimensions between Halfspace matrix and Halfspace vector');
            end

            if mH ~= obj.dim
                error('Inconsistent dimensions between Halfspace and star set');
            end

            m = size(obj.V, 2);
            C1 = H * obj.V(:, 2:m);
            d1 = g - H * obj.V(:, 1);

            new_C = vertcat(obj.C, C1);
            new_d = vertcat(obj.d, d1);

            O = Octatope(obj.V, new_C, new_d, obj.predicate_lb, obj.predicate_ub, obj.state_lb, obj.state_ub);

            if O.isEmptySet
                O = [];
            end
        end

    end
    

    methods % get methods (also estimate)

        % find range of a state at specific position
        function [xmin, xmax] = getRange(varargin)
            % @index: position of the state
            % @lp_solver: (optional) the name of the desired LP solver to use
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 11/16/2018

            switch nargin
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end
            
            f = obj.V(index, 2:obj.nVar + 1);

            if all(f(:) == 0)
                xmin = obj.V(index,1);
                xmax = obj.V(index,1);
            else
                [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                if ismember(exitflag, ["l1", "g5"])
                    xmin = fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = ' + string(exitflag));
                end
                
                [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub); 
                if ismember(exitflag, ["l1","g5"])
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution');
                end
            end
        end

        function xmin = getMin(varargin)
            % @index: position of the state
            % @lp_solver: (optional) the name of the desired LP solver to use
            % xmin: min value of x[index]

            switch nargin
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end

            if index < 1 || index > obj.dim
                error('Invalid index');
            end

            f = obj.V(index, 2:obj.nVar + 1);
            if all(f(:) == 0)
                xmin = obj.V(index, 1)
            else
                [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, lp_solver);
                if ismember(exitflag, ["l1", "g5"])
                    xmin = fval + obj.V(index, 1)
                else
                    error('Cannot find an optimal solution, exitflag = ' + exitflag);
                end
            end
        end

        function xmin = getMins(varargin)
            % @map: an array of indexes
            % xmin: min values of x[indexes]

            switch nargin
                case 5
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver  = varargin{5};
                case 4
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = [];
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    map = varargin{2}; 
                    par_option = 'single';
                    dis_option = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of inputs, should be 1, 2, 3, or 4');
            end

            n = length(map);
            xmin = cast(zeros(n, 1), 'like', obj.V);
            if isempty(par_option) || strcmp(par_option, 'single') % get mins using single core
                reverseStr = '';
                for i = 1:n
                    xmin(i) = obj.getMin(map(i), lp_solver);
                    if strcmp(dis_option, 'display')
                        msg = sprintf('%d/%d', i, n);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            elseif strcmp(par_option, 'parallel') % get mins using multiple cores
                f = obj.V(map, 2:obj.nVar + 1);
                V1 = ojb.V(map, 1);
                C1 = obj.C;
                d1 = obj.d;
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
                parfor i=1:n
                    if all(f(i, :) == 0)
                        xmin(i) = V1(i, 1);
                    else
                        [fval, exitflag] = lpsolver(f(i, :), C1, d1, [], [], pred_lb, pred_ub, lp_solver);
                        if ismember(exitflag, ["l1", "g5"])
                            xmin(i) = fval + V1(i, 1);
                        else
                            error('Cannot find an optimal solution, exitflag = %d', exitflag);
                        end
                    end
                end
            else
                error('Unknown parallel option');
            end
        end

        % get max
        function xmax = getMax(varargin)
            % @index: position of the state
            % @lp_solver: (optional) the name of the desired LP solver to use
            % xmax: max value of x[index]

            switch nargin
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end

            if index < 1 || index > obj.dim
                error('Invalid index');
            end

            f = obj.V(index, 2:obj.nVar + 1);
            if all(f(:) == 0)
                xmax = obj.V(index, 1);
            else
                [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, lp_solver);
                if ismember(exitflag, ["l1", "g5"])
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = %s', exitflag);
                end
            end
        end

        % get maxs
        function xmax = getMaxs(varargin)
            % @map: an array of indexes
            % xmax: max values of x[indexes]
            
            switch nargin
                case 5
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver  = varargin{5};
                case 4
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = [];
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    map = varargin{2}; 
                    par_option = 'single';
                    dis_option = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of inputs, should be 1, 2, 3, or 4');
            end
            
            n = length(map);
            xmax = cast(zeros(n, 1), 'like', obj.V);
                        
            if isempty(par_option) || strcmp(par_option, 'single') % get maxs using single core
                reverseStr = '';
                for i = 1:n
                    xmax(i) = obj.getMax(map(i), lp_solver);
                    if strcmp(dis_option, 'display')
                        msg = sprintf('%d/%d', i, n);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            elseif strcmp(par_option, 'parallel') % get maxs using multiple cores 
                f = obj.V(map, 2:obj.nVar + 1);
                V1 = obj.V(map, 1);
                C1 = obj.C;
                d1 = obj.d;
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
                parfor i=1:n
                    if all(f(i,:)==0)
                        xmax(i) = V1(i,1);
                    else
                        [fval, exitflag] = lpsolver(-f(i, :), C1, d1, [], [], pred_lb, pred_ub, lp_solver); 
                        if ismember(exitflag, ["l1","g5"])
                            xmax(i) = -fval + V1(i, 1);
                        else
                            error("Cannot find an optimal solution, exitflag = " + string(exitflag));
                        end                                   
                    end
                end
            else
                error('Unknown parallel option');
            end
     
        end

        % get lower bound and upper bound vector of the state variables
        function [lb, ub] = getRanges(obj)
            
            if ~obj.isEmptySet
                n = obj.dim;
                lb = zeros(n, 1);
                ub = zeros(n, 1);
                for i = 1:n
                    [lb(i), ub(i)] = obj.getRange(i);
                end
            else
                lb = [];
                ub = [];
            end
        
        end

    end


    methods % check methods

        % check if empty set
        function bool = isEmptySet(obj)
            f = zeros(1, obj.nVar, 'like', obj.V) % objective function
            [~, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, 'linprog', 'emptySet');
            if ismember(exitflag, ["l1", "g2", "g5"])
                bool = 0;
            elseif ismember(exitflag, ["l-2", "l-5", "g3", "g4", "g110"])
                bool = 1;
            else
                error('Error, exitflag = %d', exitflag);
            end
        
        end
        
    end

end

