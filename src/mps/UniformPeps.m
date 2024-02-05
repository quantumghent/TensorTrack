classdef UniformPeps
    % UniformPeps - Implementation of infinite translation invariant PEPS
    %
    
    properties
        A cell
    end
    
    
    %% Constructors
    methods
        function peps = UniformPeps(varargin)
            % Usage
            % -----
            % :code:`peps = UniformMps(A)`
            %
            % Arguments
            % ---------
            % A : :class:`cell` of :class:`PepsTensor`
            %
            % Returns
            % -------
            % peps : :class:`UniformPeps`
            
            if nargin == 0, return; end % default empty constructor
            
            if nargin == 1
                if isa(varargin{1}, 'UniformPeps') % copy constructor
                    peps.A = varargin{1}.A;

                elseif isa(varargin{1}, 'Tensor')
                    peps.A{1} = PepsTensor(varargin{1});
                    
                elseif isa(varargin{1}, 'PepsTensor')
                    peps.A{1} = varargin{1};
                    
                elseif iscell(varargin{1})
                    assert(isa(varargin{1}{1},'PepsTensor'))
                    for i = height(varargin{1}):-1:1
                        for j = width(varargin{1}):-1:1
                            peps.A{i,j} = varargin{1}{i,j};
                        end
                    end
                else
                    error('Invalid constructor for UniformPeps.')
                end
                
            else
                error('Invalid constructor for UniformPeps.')
            end
        end
    end
    
    %% Properties
    methods
        function h = height(peps)
            % vertical period over which the peps is translation invariant.
            h = height(peps.A);
        end
        
        function w = width(peps)
            % horizontal period over which the peps is translation invariant.
            w = width(peps.A);
        end
        
        function s = westvspace(peps, h, w)
            % return the virtual space to the left of site (h,w).
            if nargin == 1, h = 1:height(peps); w = 1:width(peps); end
            s = cellfun(@westvspace, peps.A(h,w));
        end

        function s = southvspace(peps, h, w)
            % return the virtual space to the bottom of site (h,w).
            if nargin == 1, h = 1:height(peps); w = 1:width(peps); end
            s = cellfun(@southvspace, peps.A(h,w));
        end

        function s = eastvspace(peps, h, w)
            % return the virtual space to the right of site (h,w).
            if nargin == 1, h = 1:height(peps); w = 1:width(peps); end
            s = cellfun(@eastvspace, peps.A(h,w));
        end

        function s = northvspace(peps, h, w)
            % return the virtual space to the top of site (h,w).
            if nargin == 1, h = 1:height(peps); w = 1:width(peps); end
            s = cellfun(@northvspace, peps.A(h,w));
        end
        
        function s = pspace(peps, h, w)
            % return the physical space at site (h,w).
            if nargin == 1, h = 1:height(peps); w = 1:width(peps); end
            s = cellfun(@pspace, peps.A(h,w));
        end

        function peps = rot90(peps)
            peps.A = cellfun(@(x)rot90(x),peps.A.','UniformOutput',false);
        end

        function peps = rot270(peps)
            peps.A = cellfun(@(x)rot270(x),peps.A.','UniformOutput',false);
        end

        function type = underlyingType(peps)
            type = underlyingType(peps.A{1,1});
        end
    end
    
    
    %% Methods
    methods
        %build horizontalTransferMatrix

        %build verticalTransferMatrix

        %CtmrgEnvironment

        function ctmrgenv = CtmrgEnvironment(peps_top, peps_bot, varargin)

            h = height(peps_top);
            w = width(peps_top);
            C = cell(4, h, w);
            T = cell(4, h, w);

            if nargin == 2
                vspace = repmat(one(space(peps_top.A{1})), h, w);
            elseif all(size(varargin{1})==[1,1])
                vspace = repmat(varargin{1}, h, w);
            else
                vspace = varargin{1};
            end
            
            for i = 1:h
                for j = 1:w
                    C{1,i,j} = Tensor.randnc(vspace(i,j), vspace(i,j));
                    C{2,i,j} = Tensor.randnc(vspace(i,prev(j,w)), vspace(i,prev(j,w)));
                    C{3,i,j} = Tensor.randnc(vspace(prev(i,h),prev(j,w)), vspace(prev(i,h),prev(j,w)));
                    C{4,i,j} = Tensor.randnc(vspace(prev(i,h),j), vspace(prev(i,h),j));

                    T{1,i,j} = Tensor.randnc([vspace(i,prev(j,w)), northvspace(peps_top,next(i,h),j)', northvspace(peps_bot,next(i,h),j)], vspace(i,j));
                    T{2,i,j} = Tensor.randnc([vspace(prev(i,h),prev(j,w)), eastvspace(peps_top,i,prev(j,w))', eastvspace(peps_bot,i,prev(j,w))], vspace(i,prev(j,w)));
                    T{3,i,j} = Tensor.randnc([vspace(prev(i,h),j), southvspace(peps_top,prev(i,h),j)', southvspace(peps_bot,prev(i,h),j)], vspace(prev(i,h),prev(j,w)));
                    T{4,i,j} = Tensor.randnc([vspace(i,j), westvspace(peps_top,i,next(j,w))', westvspace(peps_bot,i,next(j,w))], vspace(prev(i,h),j));
                end
            end

            ctmrgenv = CtmrgEnvironment(C, T);

        end

    end
end

