classdef GlobalAveragePooling3DLayer < handle
    % GlobalAveragePooling 3D Layer object
    % downsamples the input by computing the mean of the height and width dimensions
    % 
    % Author: Samuel Sasaki
    % Date: 08/27/2025
    
    properties
        Name = 'GlobalAvgPoolingLayer';          % default
        NumInputs = 1;              % default
        InputNames = {'in1'}; % default
        NumOutputs = 1;             % default
        OutputNames = {'out'};      % default
    end
    
    methods % constructor
        
        % create layer
        function obj = GlobalAveragePooling3DLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            
            switch nargin
                
                case 5
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                case 1
                    name = varargin{1};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in1'};
                    outputNames = {'out'};
                case 0
                    name = 'global_average_pooling_3d';
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in1'};
                    outputNames = {'out'};
                otherwise
                    error('Invalid number of input arguments, should be 0, 1 or 5');        
            end
            
            if ~ischar(name)
                error('Invalid name, should be a charracter array');
            end
            
            if numInputs < 1
                error('Invalid number of inputs');
            end
                       
            if numOutputs < 1
                error('Invalid number of outputs');
            end
            
            if ~iscell(inputNames)
                error('Invalid input names, should be a cell');
            end
            
            if ~iscell(outputNames)
                error('num of inputs do not match with num of input names');
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames; 
        end

        % change params to gpuArrays
        function obj = toGPU(obj)
            % nothing to change in here (no params)
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, ~)
            
        end
            
    end
        
        
    methods % main methods
        
        % evaluate
        function output = evaluate(obj, input)
            % addition layer takes usually two inputs, but allows many (N)
            %
            input = dlarray(input, "SSSCB");
            output = extractdata(avgpool(input,'global'));
        end
 
        % reach (TODO)
        function output = reach_single_input(obj, input)
            % @in_image: input volumestar
            % @image: output set
            
            if ~isa(input, 'VolumeStar')
                error('The input is not an VolumeStar object');
            end
            Y = obj.evaluate(input.V);                       
            output = VolumeStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);
            
        end
        
        % handle multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output VolumeStar
            
            n = length(inputs);
            if isa(inputs(1), 'VolumeStar')
                S(n) = VolumeStar;
            else
                error('Unknown input data set');
            end
          
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        % reachability analysis with multiple inputs
        function VS = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
           
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = varargin{3};
                    option = varargin{4}; % computation option

                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; % don't care the rest inputs
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5 or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || strcmp(method, 'approx-zono') || contains(method, "relax-star")
                VS = obj.reach_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
  
        end
        
        function VS = reachSequence(varargin)
            obj = varargin{1};
            VS = obj.reach(varargin{2:end});
        end
    end
    
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer from matlab
                      
            if ~isa(layer, 'nnet.cnn.layer.GlobalAveragePooling3DLayer') 
                error('Input is not a GlobalAveragePooling3DLayer layer');
            end
            L = GlobalAveragePooling3DLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames);
        end

    end
end
