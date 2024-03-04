function [sequence, cost] = netcon(legLinks,verbosity,costType,muCap,allowOPs,legCosts)
% Finds most efficient way to contract a tensor network
% v2.00 by Robert N. C. Pfeifer 2014
% Contact: rpfeifer.public@gmail.com, robert.pfeifer@mq.edu.au
%
% Usage
% -----
% :code:`[sequence, cost] = netcon(legLinks, verbosity, costType, muCap, allowOPs, legCosts)`
%
% Arguments
% ---------
% legLinks
%   abelling of tensor indices.
%
% Optional Arguments
% ------------------
% verbosity
%   0: Quiet. 1: State final result. 2: State intermediate and final results. 3: Also display object sizes during pairwise contractions. (2)
% costType
%   1: Absolute value. 2: Multiple and power of chi. (2)
% muCap
%   Initially restrict search to sequences having a cost of muCap (if costType==1) or
%   O(X^muCap) (if costType==2). This value will increment automatically if required, and it
%   is recommended that it be left at the default value of 1. (1)
% allowOPs
%   Allow contraction sequences including outer products: 0/false: No. 1/true: Yes. (true)
% legCosts
%   For costType==1: nx2 table. A row reading [a b] assigns a dimension of b to index a.
%   Default: 2 for all legs. For costType==2: nx3 table. A row reading [a b c] assigns a
%   dimension of bX^c to index a, for unspecified X. Default: 1X^1 for all legs.

arguments
    legLinks
    verbosity   = 2 % 0: No displayed output. 1: Display final result. 2: Display progress reports.
    costType    = 2 % 1: Absolute value. 2: Multiple and power of chi.
    muCap       = 1
    allowOPs    = 1
    legCosts    = []
end

% Benchmarking examples
% ---------------------
% 3:1 1D MERA:
% tic;netcon({[-1 1 2 3],[2 4 5 6],[1 5 7 -3],[3 8 4 9],[6 9 7 10],[-2 8 11 12],[10 11 12 -4]},0,2,1,1);toc
% All in MATLAB, no OPs: 0.041s
% All in MATLAB, with OPs: 0.054s
% Using C++, no OPs: 0.0019s
% Using C++, with OPs: 0.0019s
%
% 9:1 2D MERA:
% tic;netcon({[1 4 5 6 7 8],[2 9 10 11 12 13],[3 14 15 16 17 18],[-6 9 23 24 25 26],[-5 5 19 20 21 22],[8 14 27 28 29 30],[12 15 31 32 33 34],[22 25 28 31 35 36 37 38],[-4 20 23 35 39 40 41 42],[42 36 37 38 43 44 45 46],[41 24 44 26 47 48],[19 40 21 43 49 50],[27 45 29 30 51 52],[46 32 33 34 53 54],[-2 -3 39 49 47 55],[4 50 6 7 51 56],[48 10 11 53 13 57],[52 54 16 17 18 58],[55 56 57 58 -1 1 2 3]},0,2,1,1);toc
% All in MATLAB, no OPs: 5.4s
% All in MATLAB, with OPs: 6.1s
% Using C++, no OPs: 0.066s
% Using C++, with OPs: 0.069s
%
% 4:1 2D MERA:
% tic;netcon({[64 72 73 19 22],[65 74 75 23 76],[66 77 20 79 28],[67 21 24 29 34],[68 25 78 35 80],[69 81 30 83 84],[70 31 36 85 86],[71 37 82 87 88],[-5 19 20 21 1 2 4 5],[22 23 24 25 3 26 6 27],[28 29 30 31 7 8 32 33],[34 35 36 37 9 38 39 40],[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18],[10 11 13 14 41 42 43 44],[12 26 15 27 47 48 49 50],[16 17 32 33 53 54 55 56],[18 38 39 40 60 61 62 63],[-2 -3 -4 41 89],[72 73 42 47 90],[74 75 48 76 45],[77 43 79 53 46],[44 49 54 60 51],[50 78 61 80 52],[81 55 83 84 57],[56 62 85 86 58],[63 82 87 88 59],[89 90 45 46 51 52 57 58 59 -1 64 65 66 67 68 69 70 71]},0,2,1,1);toc
% All in MATLAB, no OPs: 2486s (41.4m)
% All in MATLAB, with OPs: 3919.6s (65.3m)
% Using C++, no OPs: 36.12s
% Using C++, with OPs: 36.49s


% Changes in v2.00:
% -----------------
% Changed algorithm to breadth-first search instead of dynamic programming.
% Made handling of subnetworks iterative, not recursive.
% Made displayed output more user-friendly for networks composed of disjoint subnetworks.
% Added warning, not error, for tensors or subnetworks of trivial dimension (i.e. just a number).
% Added more vigorous exclusion of unnecessary outer products.

% Changes in v1.01:
% -----------------
% Improved structure of code.
% Fixed bug returning wrong cost when trivial indices present and non-consecutive numbering of positive indices.
% Removed keepTriv option.
% Fixed omission of trailing zeros in sequence when network is disjoint and all contractions are traces.
% Corrected helptext.


% Check input data (and strip trivial indices from legLinks)
% ================
if ~isempty(legCosts)
    [trivialindices legCosts legLinks] = checkInputData(legLinks,verbosity,costType,muCap,allowOPs,legCosts);
else
    [trivialindices legCosts legLinks] = checkInputData(legLinks,verbosity,costType,muCap,allowOPs);
end

% Divide network into disjoint subnets
% ====================================
subnetlist = zeros(1,numel(legLinks));
subnetcounter = 1;
while any(subnetlist==0)
    flags = findSubnetIncluding(legLinks,find(subnetlist==0,1,'first'));
    subnetlist(flags) = subnetcounter;
    subnetcounter = subnetcounter + 1;
end
subnetcounter = subnetcounter-1;

% Evaluate contraction sequences for disjoint subnets
% ===================================================
sequence = cell(1,subnetcounter);
cost = cell(1,subnetcounter);
freelegs = cell(1,subnetcounter);
donetrace = false(1,subnetcounter);
for a=1:subnetcounter
    if verbosity > 0
        disp(' ');
        if subnetcounter > 1
            disp(['Subnet ' num2str(a) ':']);
        end
        if verbosity > 1
            disp(['Tensors: ' unpaddednum2str(find(subnetlist==a))]);
        end
    end
    [sequence{a} cost{a} freelegs{a} donetrace(a)] = netcon_getsubnetcost(legLinks(subnetlist==a),verbosity,costType,muCap,legCosts,allowOPs);
end
donetrace = any(donetrace);

% Perform outer products of subnets, add costs, and merge sequences
% =================================================================
[sequence cost] = performOPs(sequence,cost,freelegs,legCosts,costType,verbosity);
sequence = [sequence trivialindices]; % Append tracing over trivial indices on final object

% Display result
% ==============
if verbosity>0
    if ~isempty(trivialindices) && verbosity > 1
        disp(' ');
        disp(['Restore summed trivial indices: ' unpaddednum2str(trivialindices)]);
    end
    disp(' ');
    if verbosity > 1 && subnetcounter > 1
        disp('Entire network:');
    end
    netcon_displayresult(sequence,cost,donetrace,costType);
end
end

function [sequence cost] = performOPs(sequence,cost,freelegs,legCosts,costType,verbosity)
% Renumber freelegs by position and trim index label list off legCosts
for a=1:numel(freelegs)
    for b=1:numel(freelegs{a})
        freelegs{a}(b) = find(legCosts(:,1)==freelegs{a}(b));
    end
end
legCosts = legCosts(:,2:end);
% Perform outer products
if costType==2
    tensormulsize = ones(size(freelegs));
    tensorchisize = zeros(size(freelegs));
    for a=1:numel(freelegs)
        for b=1:numel(freelegs{a})
            tensormulsize(a) = tensormulsize(a)*legCosts(freelegs{a}(b),1);
            tensorchisize(a) = tensorchisize(a)+legCosts(freelegs{a}(b),2);
        end
    end
    if any(tensormulsize==1 & tensorchisize==0) && numel(freelegs)>1
        if sum(tensormulsize==1 & tensorchisize==0)==1
            if verbosity > 0
                disp(' ');
            end
            warning('netcon:trivialsubnet',['Subnet ' unpaddednum2str(find(tensormulsize==1 & tensorchisize==0)) ' reduces to a single number. Since this subnet is not the entire network, contraction sequence may not be optimal.']);
        else
            if verbosity > 0
                disp(' ');
            end
            warning('netcon:trivialsubnet',['More than one subnet reduces to a single number: Contraction sequence may not be optimal. List of subnets reducing to a single number: ' unpaddednum2str(find(tensormulsize==1 & tensorchisize==0))]);
        end
    end
    while numel(tensorchisize)>1
        % Sort tensors
        ptr = 1;
        while (ptr<numel(tensorchisize))
            if tensorchisize(ptr)>tensorchisize(ptr+1) || (tensorchisize(ptr)==tensorchisize(ptr+1) && tensormulsize(ptr)>tensormulsize(ptr+1))
                t = tensorchisize(ptr+1); tensorchisize(ptr+1) = tensorchisize(ptr); tensorchisize(ptr) = t;
                t = tensormulsize(ptr+1); tensormulsize(ptr+1) = tensormulsize(ptr); tensormulsize(ptr) = t;
                t = sequence{ptr}; sequence{ptr} = sequence{ptr+1}; sequence{ptr+1} = t;
                t = cost{ptr}; cost{ptr} = cost{ptr+1}; cost{ptr+1} = t;
                ptr = ptr - 1;
                if ptr==0
                    ptr = 2;
                end
            else
                ptr = ptr + 1;
            end
        end
        % Combine smallest two tensors
        outerproductmul = tensormulsize(1) * tensormulsize(2);
        outerproductpower = tensorchisize(1) + tensorchisize(2);
        if numel(cost{1})<outerproductpower+1
            cost{1}(outerproductpower+1) = outerproductmul;
        else
            cost{1}(outerproductpower+1) = cost{1}(outerproductpower+1) + outerproductmul;
        end
        tensormulsize = [outerproductmul tensormulsize(3:end)];
        tensorchisize = [outerproductpower tensorchisize(3:end)];
        cost{1} = addCosts(cost{1},cost{2},costType);
        sequence{1} = [sequence{1} sequence{2}];
        cost = cost([1 3:end]);
        sequence = sequence([1 3:end]);
    end
else
    tensorsize = ones(size(freelegs));
    for a=1:numel(freelegs)
        for b=1:numel(freelegs{a})
            tensorsize(a) = tensorsize(a)*legCosts(freelegs{a}(b),1);
        end
    end
    while numel(tensorsize)>1
        % Sort tensors
        tensorsize = sort(tensorsize,'ascend');
        % Combine smallest two tensors
        outerproduct = tensorsize(1) * tensorsize(2);
        cost{1} = cost{1} + outerproduct;
        tensorsize = [outerproduct tensorsize(3:end)];
        cost{1} = cost{1} + cost{2};
        sequence{1} = [sequence{1} sequence{2}];
        cost = cost([1 3:end]);
        sequence = sequence([1 3:end]);
    end
end
sequence = sequence{1};
sequence = [sequence zeros(1,numel(freelegs)-1)]; % Add outer products between disjoint subnets to the contraction sequence
cost = cost{1};
end

function [sequence cost negindices donetrace] = netcon_getsubnetcost(legLinks,verbosity,costType,muCap,legCosts,allowOPs)
% Performs any traces on single tensors
% Checks if network is trivial
% If the network is not trivial, invokes netcon_nondisj to evaluate cost for the network

% Get index lists for this subnet and trim legCosts
% =================================================
[posindices negindices] = getIndexLists(cell2mat(legLinks));

% Any contraction to do?
% ======================
if numel(legLinks)<2
    % List of tensors has only one entry
    % ==================================
    if ~isempty(posindices)
        % One tensor, with traces - do them, then go to "Finished contracting"
        if costType==2
            cost = [];
        else
            cost = 0;
        end
        sequence = posindices;
        donetrace = true;
        if verbosity > 0
            disp('Tracing costs only');
        end
        % Go to "Finished contracting"
    else
        % Nothing to do
        sequence = [];
        if costType==2
            cost = [];
        else
            cost = 0;
        end
        donetrace = false;
        if verbosity > 0
            disp('Network is trivial');
        end
        % One tensor, no traces - go to "Finished contracting"
    end
else
    % >1 tensor, tracing and/or contraction to do
    % ===========================================
    % Make positive indices consecutive
    % Make negative indices consecutive following on from positive indices
    % Truncate cost table and remove indexing column
    % (Note original labelling for legs may be recovered from posindices and negindices)
    [legLinks legCosts] = reprocessLegLinks(legLinks,posindices,negindices,legCosts);
    
    % Determine any traces (legs with both ends on a single tensor) and eliminate
    % ===========================================================================
    [tracedindices donetrace legLinks legCosts linkposindices] = eliminateTraces(legLinks,posindices,legCosts);
    
    if isempty(linkposindices)
        % No positive indices after tracing
        sequence = posindices(tracedindices);
        if costType==2
            cost = [];
        else
            cost = 0;
        end
    else
        % Find best sequence
        % ==================
        % Invoke sequence finder for traceless non-disjoint network
        for a=1:numel(legLinks)
            legLinks{a} = int32(abs(legLinks{a}));
        end
        legCosts = double(legCosts);
        verbosity = double(verbosity);
        costType = double(costType);
        muCap = double(muCap);
        allowOPs = double(allowOPs);
        posindices = double(posindices);
        tracedindices = int32(tracedindices);
        [sequence cost] = netcon_nondisj(legLinks,legCosts,verbosity,costType,muCap,allowOPs,posindices,tracedindices);
%         try
%             [sequence cost] = netcon_nondisj_cpp(legLinks,legCosts,verbosity,costType,muCap,allowOPs,posindices,tracedindices);
%         catch ME
%             if isequal(ME.identifier,'MATLAB:UndefinedFunction')
%                 [sequence cost] = netcon_nondisj(legLinks,legCosts,verbosity,costType,muCap,allowOPs,posindices,tracedindices);
%             else
%                 rethrow(ME);
%             end
%         end
        % Reinsert indices traced at the beginning of the contraction process
        sequence = [tracedindices sequence];
        % Translate sequence back into user's original leg labels
        sequence(sequence>0) = posindices(sequence(sequence>0));
    end
end
end

function [tracedindices donetrace legLinks legCosts posindices] = eliminateTraces(legLinks,posindices,legCosts)
tracedindices = [];
donetrace = false;
% For each tensor:
for a=1:numel(legLinks)
    % - Find all tracing legs on that tensor, note the indices, and delete them.
    % - Tracing can be performed over all traced indices on a tensor simultaneously by combining them - no need to sort
    b = 1;
    while b<numel(legLinks{a})
        if any(legLinks{a}(b+1:end)==legLinks{a}(b))
            tracedindices(end+1) = legLinks{a}(b); %#ok<AGROW>
            legLinks{a}(legLinks{a}==legLinks{a}(b)) = [];
        else
            b = b + 1;
        end
    end
end
% Then:
if ~isempty(tracedindices)
    donetrace = true;
    % - Delete removed rows from leg cost list
    legCosts(tracedindices,:) = [];
    % - Renumber all positive and negative legs consecutively again
    renumindices = 1:numel(posindices);
    renumindices(tracedindices) = [];
    posindices(tracedindices) = [];
    for a=1:numel(legLinks)
        for b=1:numel(legLinks{a})
            if legLinks{a}(b)>0
                legLinks{a}(b) = find(renumindices==legLinks{a}(b));
            else
                legLinks{a}(b) = legLinks{a}(b) + numel(tracedindices);
            end
        end
    end
end
end

function [legLinks legCosts] = reprocessLegLinks(legLinks,posindices,negindices,legCosts)
% Renumber all positive and negative indices consecutively
% ========================================================
for a=1:numel(legLinks)
    for b=find(legLinks{a}>0)
        legLinks{a}(b) = find(posindices==legLinks{a}(b));
    end
    for b=find(legLinks{a}<0)
        legLinks{a}(b) = -find(negindices==legLinks{a}(b))-numel(posindices);
    end
end

% Assemble truncated cost table
% =============================
for a=size(legCosts,1):-1:1
    if ~any([posindices negindices]==legCosts(a,1))
        legCosts(a,:) = [];
    end
end
legCosts = legCosts(:,2:end);
end

function cost = addCosts(cost1,cost2,costType)
if costType==1
    cost = cost1 + cost2;
else
    cost = zeros(1,max([numel(cost1) numel(cost2)]));
    cost(1:numel(cost1)) = reshape(cost1,1,numel(cost1));
    cost(1:numel(cost2)) = cost(1:numel(cost2)) + reshape(cost2,1,numel(cost2));
end
end

function flags = findSubnetIncluding(legLinks,pos)
new = true;
flags = zeros(1,numel(legLinks)); % 0: May not be in subnet
flags(pos) = 1; % 1: Follow connections off this tensor
while new
    new = false;
    fromlist = find(flags==1);
    flags(flags==1) = 2; % 2: Have followed all connections from this tensor
    for checkfrom = fromlist
        for checkto = find(flags==0)
            for x = 1:numel(legLinks{checkfrom})
                if any(legLinks{checkto}==legLinks{checkfrom}(x))
                    flags(checkto) = true;
                    new = true;
                end
            end
        end
    end
end
flags = flags~=0;
end

function [trivialindices legCosts legLinks] = checkInputData(legLinks,verbosity,costType,muCap,allowOPs,legCosts)
% Check data sizes
if ~iscell(legLinks)
    error('legLinks must be a cell array')
end
if numel(legLinks)==0
    error('legLinks may not be empty')
end
if ~isnumeric(verbosity) || ~any([0 1 2 3]==verbosity)
    error('verbosity must be 0, 1, 2 or 3')
end
if ~isnumeric(muCap) || numel(muCap)~=1 || ~isreal(muCap) || muCap<=0
    error('muCap must be a real positive number');
end
if ~isnumeric(costType) || ~any([1 2]==costType)
    error('costType must be either 1 or 2')
end
if (~isnumeric(allowOPs) && ~islogical(allowOPs)) || ~any([0 1]==allowOPs)
    error('allowOPs must be either 0 or 1');
end
for a=1:numel(legLinks)
    if ~isnumeric(legLinks{a})
        error('Entries in legLinks must be numeric arrays')
    end
end
if size(legLinks,1)~=1 || size(legLinks,2)~=numel(legLinks)
    error('Array of index connections (legLinks) has incorrect dimension - should be 1xn')
end
for a=1:numel(legLinks)
    if size(legLinks{a},1)~=1 || size(legLinks{a},2)~=numel(legLinks{a})
        error(['legLinks entry ' num2str(a) ' has wrong dimension - should be 1xn']);
    end
    if isempty(legLinks{a})
        error(['Empty list of indices on tensor ' num2str(a)])
    end
end
allindices = cell2mat(legLinks);
if any(allindices==0)
    error('Zero entry in index list')
elseif any(imag(allindices)~=0)
    error('Complex entry in index list')
elseif any(int32(allindices)~=allindices)
    error('Non-integer entry in legLinks');
end
posindices = sort(allindices(allindices>0),'ascend');
negindices = sort(allindices(allindices<0),'descend');
if ~isempty(posindices)
    % Test all positive indices occur exactly twice
    if mod(numel(posindices),2)~=0
        maxposindex = posindices(end);
        posindices = posindices(1:end-1);
    end
    flags = (posindices(1:2:numel(posindices))-posindices(2:2:numel(posindices)))~=0;
    if any(flags)
        errorpos = 2*find(flags~=0,1,'first')-1;
        if errorpos>1 && posindices(errorpos-1)==posindices(errorpos)
            error(['Error in index list: Index ' num2str(posindices(errorpos)) ' appears more than twice']);
        else
            error(['Error in index list: Index ' num2str(posindices(errorpos)) ' only appears once']);
        end
    end
    if exist('maxposindex','var')
        if posindices(end)==maxposindex
            error(['Error in index list: Index ' num2str(maxposindex) ' appears more than twice']);
        else
            error(['Error in index list: Index ' num2str(maxposindex) ' only appears once']);
        end
    end
    posindices = posindices(1:2:numel(posindices)); % List of all positive indices, sorted ascending, each appearing once.
    flags = posindices(1:end-1)==posindices(2:end);
    if any(flags)
        errorpos = find(flags,1,'first');
        error(['Error in index list: Index ' num2str(posindices(errorpos)) ' appears more than twice']);
    end
else
    posindices = [];
end
% Test all negative indices occur exactly once
flags = negindices(1:end-1)==negindices(2:end);
if any(flags)
    errorpos = find(flags,1,'first');
    error(['Error in index list: Index ' num2str(negindices(errorpos)) ' appears more than once']);
end

% Check leg size data
trivialindices = [];
if ~exist('legCosts','var')
    if costType==1
        % Create basic leg costs list (all 2)
        legCosts = [[posindices.';negindices.'] 2*ones(numel(posindices)+numel(negindices),1)];
    else
        % Create basic leg costs list (all Chi^1)
        legCosts = [[posindices.';negindices.'] ones(numel(posindices)+numel(negindices),2)];
    end
else
    if ~isnumeric(legCosts)
        error('legCosts must be numeric')
    end
    % Check valid leg costs list, & process
    if ndims(legCosts)>2 || any(legCosts(:,1)==0) || size(legCosts,2)~=2+(costType==2) || any(legCosts(:,2)<=0)
        error('Bad index dimensions list (error in legCosts: type ''help netcon'' for specifications)')
    end
    if costType==2
        if any(legCosts(:,3)<0)
            error('For index dimensions of the form aX^b, values of b less than 0 are not supported')
        end
        if any(legCosts(:,2)<=0)
            error('For index dimensions of the form aX^b, values of a less than or equal to 0 are not supported')
        end
    else
        if any(legCosts(:,2)<=0)
            error('Index dimensions less than or equal to 0 are not supported')
        end
    end
    [t1 ix1] = sort(legCosts(legCosts(:,1)>0,1));
    [t2 ix2] = sort(legCosts(legCosts(:,1)<0,1),'descend');
    t1 = reshape(t1,1,[]);
    t2 = reshape(t2,1,[]);
    flag = t1(1:end-1)==t1(2:end);
    if any(flag)
        error(['Dimension of index ' num2str(t1(find(flag,1))) ' specified twice'])
    end
    if ~(isempty(t1) && isempty(posindices)) % Necessary as a 1x0 matrix is not the same as a 0x0 matrix
        if ~isequal(t1,posindices)
            error(['Dimension not specified for all positive indices. Supplied: ' num2str(t1) '     Required: ' num2str(posindices)])
        end
    end
    flag = t2(1:end-1)==t2(2:end);
    if any(flag)
        error(['Dimension of index ' num2str(t2(find(flag,1))) ' specified twice'])
    end
    if ~(isempty(t2) && isempty(negindices))
        if ~isequal(t2,negindices)
            error(['Dimension not specified for all negative indices. Supplied: ' num2str(t2) '     Required: ' num2str(negindices)])
        end
    end
    
    % Store list of any summed trivial indices to be stripped
    if costType==1
        trivialindices = legCosts(legCosts(:,2)==1,1);
    else
        trivialindices = legCosts(legCosts(:,2)==1 & legCosts(:,3)==0,1);
    end
    trivialindices = reshape(sort(trivialindices(trivialindices>0),'descend'),1,[]);
    
    % Strip trivial indices
    for a=numel(trivialindices):-1:1
        for b=1:numel(legLinks)
            if sum(legLinks{b}==trivialindices(a))==2
                trivialindices(a) = []; % Do not strip trivial single-tensor traces (not that this actually matters either way)
            else
                legLinks{b}(legLinks{b}==trivialindices(a)) = [];
            end
        end
    end
    if ~isempty(trivialindices) && verbosity > 1
        disp(' ');
        disp(['Ignore summed trivial indices: ' unpaddednum2str(trivialindices)]);
    end
    
    % Order leg costs list
    t1 = legCosts(legCosts(:,1)>0,2:end);
    t1 = t1(ix1,:);
    t2 = legCosts(legCosts(:,1)<0,2:end);
    t2 = t2(ix2,:);
    legCosts = [[posindices.';negindices.'] [t1;t2]];
end
end

function [posindices negindices] = getIndexLists(allindices)
posindices = sort(allindices(allindices>0),'ascend');
negindices = sort(allindices(allindices<0),'descend');
posindices = posindices(1:2:end);
end

function netcon_displayresult(sequence,cost,donetrace,costType)
% Displays nicely-formatted final output
t = ['Best sequence:  ' unpaddednum2str(sequence)];
if isempty(sequence)
    t = [t '<empty>'];
end
disp(t);
if isempty(cost) || isequal(cost,0)
    if donetrace
        disp('Cost:           Tracing costs only');
    else
        disp('Cost:           0');
    end
else
    t = 'Cost:           ';
    if costType==2
        for a=numel(cost):-1:2
            t = [t num2str(cost(a)) 'X^' num2str(a-1) ' + ']; %#ok<AGROW>
        end
        t = [t num2str(cost(1)) 'X^0'];
    else
        t = [t num2str(cost)];
    end
    if donetrace
        t = [t ' + tracing costs'];
    end
    disp(t);
end
end

function str = unpaddednum2str(row)
str = [];
for a=row(1:end-1)
    str = [str num2str(a) ' ']; %#ok<AGROW>
end
if ~isempty(row)
    str = [str num2str(row(end))];
end
end

% Functions used in the pure-MATLAB version only:

function [sequence cost] = netcon_nondisj(legLinks,legCosts,verbosity,costType,muCap,allowOPs,posindices,tracedindices)
% Set up initial data
% ===================
numtensors = numel(legLinks);
if costType==1
    zerocost = 0;
else
    zerocost = [];
end
allowOPs = (allowOPs==1);

% This code uses the nomenclature of Appendix G, in which a tensor which may be contracted with the result of an outer product [e.g. C in Fig.4(a)] is denoted X.

% Structure of "objects": objects{numElements}{positionInList}{legFlags,tensorFlags,sequenceToBuild,costToBuild,isOP,OPmaxdim,allIn}
% legFlags: Legs present on object
% tensorFlags: Tensor constituents of object
% sequenceToBuild: Cheapest identified contraction sequence for constructing this object
% costToBuild: Minimum identified cost of constructing this object
% isOP: Indicates whether the contraction which yielded this object was an outer product
% OPmaxdim: If isOP==true, then OPmaxdim gives the dimension of the largest constituent tensor
% allIn: If this tensor is an object X which may be contracted with an outer product, allIn is the dimension of the tensor which contributed no external legs to X, i.e. xi_C in Fig.5(c).

objects = cell(1,numtensors);
objects{1} = cell(1,numtensors);
tensorflags = false(1,numtensors);
numleglabels = size(legCosts,1);
legflags = false(1,numleglabels);

% ### Create lists used in enforcing Sec. II.B.2.c (structure of tensor X contractable with an outer product)
tensorXlegs = {};
tensorXflags = [];
tensorXdims = {};
% ### End of creating lists

for a=1:numtensors
    tensorflags(a) = true;
    legflags(abs(legLinks{a})) = true;
    objects{1}{a} = {legflags,tensorflags,[],zerocost,false,[],zeros(1,costType)};
    
    % ### Set up initial list data used in enforcing Sec. II.B.2.c
    if allowOPs
        [tensorXlegs tensorXdims tensorXflags] = addToTensorXlist(tensorXlegs,tensorXdims,tensorXflags,legflags,Inf*ones(1,costType),costType,true);
    end
    % ### End of setting up linked list data
    
    tensorflags(a) = false;
    legflags(abs(legLinks{a})) = false;
end

% ### Set up initial list data used in enforcing Sec. II.B.2.c (continued)
tensorXflags = zeros(1,numel(tensorXflags));
% ### End of setting up initial linked list data (continued)

for a=2:numtensors
    objects{a} = {};
end

newobjectflags = cell(1,numtensors);
newobjectflags{1} = true(1,numtensors);

oldmuCap = 0;
newmuCap = Inf;

done = false;
while ~done
    if verbosity>0
        if oldmuCap~=muCap
            if costType==1
                disp(['Looking for solutions with maximum cost of ' num2str(muCap)]);
            else
                disp(['Looking for solutions of cost O(X^' num2str(muCap) ')']);
            end
        end
    end
    for numInObjects = 2:numtensors
        if verbosity > 2
            disp(['Pairwise contractions (AB) involving ' num2str(numInObjects) ' fundamental tensors:']);
        end
        for numInPieceOne = 1:floor(numInObjects/2)
            numInPieceTwo = numInObjects - numInPieceOne;
            if verbosity > 2
                disp(['A contains ' num2str(numInPieceOne) ', B contains ' num2str(numInPieceTwo)'.']);
                pause(0.01);
            end
            if numel(objects{numInPieceOne})>0 && numel(objects{numInPieceTwo})>0
                
                % Iterate over pairings: Iterate over object 1
                for a = 1:numel(objects{numInPieceOne})
                    
                    % Get data for object 1
                    obj1data = objects{numInPieceOne}{a};
                    % obj1data: legs, tensors, seqToBuild, costToBuild, isOP, maxdimOP, allIn
                    legs1 = obj1data{1};
                    tensorsIn1 = obj1data{2};
                    seq1 = obj1data{3};
                    costToBuild1 = obj1data{4};
                    isOP1 = obj1data{5};
                    OPmaxdim1 = obj1data{6};
                    allIn1 = obj1data{7};
                    isnew1 = newobjectflags{numInPieceOne}(a);
                    
                    % Iterate over pairings: Iterate over object 2
                    if numInPieceOne == numInPieceTwo
                        obj2list = a+1:numel(objects{numInPieceTwo});
                    else
                        obj2list = 1:numel(objects{numInPieceTwo});
                    end
                    for b = obj2list
                        
                        % Check object 1 and object 2 don't share any common tensors (which would then appear twice in the resulting network)
                        if ~any(objects{numInPieceOne}{a}{2} & objects{numInPieceTwo}{b}{2})
                            
                            % Get data for object 2
                            obj2data = objects{numInPieceTwo}{b};
                            % obj2data: legs, tensors, seqToBuild, costToBuild, isOP, maxdimOP, allIn
                            legs2 = obj2data{1};
                            tensorsIn2 = obj2data{2};
                            seq2 = obj2data{3};
                            costToBuild2 = obj2data{4};
                            isOP2 = obj2data{5};
                            OPmaxdim2 = obj2data{6};
                            allIn2 = obj2data{7};
                            isnew2 = newobjectflags{numInPieceTwo}(b);
                            
                            commonlegs = legs1 & legs2; % Common legs
                            freelegs = xor(legs1,legs2);
                            freelegs1 = legs1 & freelegs;
                            freelegs2 = legs2 & freelegs;
                            commonlegsflag = any(commonlegs);
                            
                            isOK = allowOPs || commonlegsflag; % Exclude outer products if allowOPs is not set
                            
                            thisTensorXflag = -1;
                            % ### Enforce Sec. II.B.2.b,c,d (only perform outer product if there is a known tensor X with appropriate structure; only contract resulting object in another outer product or with an appropriate tensor X; enforce index dimension constraints)
                            if isOK && ~commonlegsflag
                                % It's an outer product. Check if a suitable tensor X exists yet to contract with this outer product [Fig.5(c) & Eq.(25)].
                                if oldmuCap == muCap
                                    % Pass for new X's
                                    if isnew1 || isnew2
                                        % Using a new object - allowed to contract with old and new X's
                                        % Note - the while/end construction is faster than using "find"
                                        Xstart = 1;
                                        Xend = 1;
                                        while Xend<=numel(tensorXflags) && tensorXflags(Xend)~=2;
                                            Xend = Xend + 1;
                                        end
                                    else
                                        % Made from old objects - only allowed to contract with new X's. Already had a chance to contract with old X's
                                        % on the previous pass so repeating these is unnecessary.
                                        Xstart = 1;
                                        while Xstart<=numel(tensorXflags) && tensorXflags(Xstart)~=1;
                                            Xstart = Xstart + 1;
                                        end
                                        Xend = Xstart;
                                        while Xend<=numel(tensorXflags) && tensorXflags(Xend)==1;
                                            Xend = Xend + 1;
                                        end
                                    end
                                else
                                    % Old X's only on this pass
                                    Xstart = 1;
                                    Xend = 1;
                                    while Xend<=numel(tensorXflags) && tensorXflags(Xend)==0;
                                        Xend = Xend + 1;
                                    end
                                end
                                for x = Xstart:Xend-1
                                    if all(tensorXlegs{x}(freelegs))
                                        % IIB2c: xi_C > xi_A (25)
                                        if isGreaterThan_sd(tensorXdims{a},getProdLegDims(freelegs,legCosts,costType),costType)
                                            % IIB2b: xi_C > xi_D && xi_C > xi_E: (16)
                                            if isGreaterThan_sd(getProdLegDims(tensorXlegs{x}&~freelegs,legCosts,costType),getProdLegDims(freelegs1,legCosts,costType),costType)
                                                if isGreaterThan_sd(getProdLegDims(tensorXlegs{x}&~freelegs,legCosts,costType),getProdLegDims(freelegs2,legCosts,costType),costType)
                                                    thisTensorXflag = tensorXflags(x);
                                                    break;
                                                end
                                            end
                                        end
                                    end
                                end
                                isOK = thisTensorXflag~=-1;
                            end
                            
                            % If either constituent is the result of an outer product, check that it is being contracted with an appropriate tensor
                            % [either this is a contraction over all indices, or this is an outer product with a tensor of larger total dimension than
                            % either constituent of the previous outer product, and Eqs. (16), (25), and (26) are satisfied].
                            if isOK && (isOP1 || isOP2)
                                % Post-OP. This contraction only allowed if it is also an outer product, or if only one object is an outer product, it involves all indices on that tensor, and the other object satisfies the relevant conditions.
                                isOK = xor(isOP1,isOP2) || ~commonlegsflag; % If contracting over common indices, only one object may be an outer product
                                if isOK
                                    if commonlegsflag
                                        % This contraction is not itself an outer product
                                        % Conditions on outer product object:
                                        if isOP1
                                            % Check all-indices condition:
                                            if isequal(commonlegs,legs1)
                                                % Check free legs on contracting tensor are larger than summing legs going to each component of outer product [Eq. (16)]
                                                isOK = isGreaterThan_sd(getProdLegDims(freelegs,legCosts,costType),OPmaxdim1,costType); % IIB2b: xi_C > xi_D, xi_C > xi_E (16)
                                            else
                                                isOK = false;
                                            end
                                        else
                                            % Check all-indices condition:
                                            if isequal(commonlegs,legs2)
                                                % Check free legs on contracting tensor are larger than summing legs going to each component of outer
                                                % product [Eq. (16)]
                                                isOK = isGreaterThan_sd(getProdLegDims(freelegs,legCosts,costType),OPmaxdim2,costType); % IIB2b: xi_C > xi_D, xi_C > xi_E (16)
                                            else
                                                isOK = false;
                                            end
                                        end
                                        % Conditions on X: Ensure X is fundamental or acceptably-constructed (note: structure is checked by requiring
                                        % non-zero value of allIn)
                                        if isOK
                                            if isOP1
                                                % Tensor 2 is X
                                                if numInPieceTwo > 1
                                                    % Tensor 2 is not fundamental
                                                    % Check tensor 2 is constructed in an acceptable fashion [Fig. 5(c) and Eqs. (25) and (26)]
                                                    isOK = isGreaterThan_sd(allIn2,getProdLegDims(freelegs,legCosts,costType),costType); % IIB2c: xi_C > xi_D (26)
                                                    if isOK
                                                        isOK = isGreaterThan_sd(allIn2,getProdLegDims(freelegs1,legCosts,costType),costType); % IIB2c: xi_C > xi_A (25)
                                                    end
                                                end
                                            else
                                                % Tensor 1 is X
                                                if numInPieceOne > 1
                                                    % Tensor 1 is not fundamental
                                                    % Check tensor 1 is constructed in an acceptable fashion [Fig. 5(c) and Eqs. (25) and (26)]
                                                    isOK = isGreaterThan_sd(allIn1,getProdLegDims(freelegs,legCosts,costType),costType); % IIB2c: xi_C > xi_D (26)
                                                    if isOK
                                                        isOK = isGreaterThan_sd(allIn1,getProdLegDims(freelegs2,legCosts,costType),costType); % IIB2c: xi_C > xi_A (25)
                                                    end
                                                end
                                            end
                                        end
                                    else
                                        % This contraction is an outer product. If either constituent is an outer product, check that both tensors
                                        % within that object are not larger than the third tensor with which they are now being contracted.
                                        if isOP1
                                            isOK = ~isGreaterThan_sd(OPmaxdim1,getProdLegDims(freelegs2,legCosts,costType),costType); % IIB2b: xi_C >= xi_A, xi_C >= xi_B (20)
                                        end
                                        if isOK && isOP2
                                            isOK = ~isGreaterThan_sd(OPmaxdim2,getProdLegDims(freelegs1,legCosts,costType),costType); % IIB2b: xi_C >= xi_A, xi_C >= xi_B (20)
                                        end
                                    end
                                end
                            end
                            % ### End of enforcing Sec. II.B.2.b,c,d  (only perform outer product if there is a known tensor X with appropriate structure; only contract resulting object in another outer product or with an appropriate tensor X; enforce index dimension constraints)
                            
                            % If contraction is not prohibited, check cost is acceptable (<=muCap and, if not involving new objects, >oldmuCap)
                            if isOK
                                % If constructing an outer product which may contract with a new X, do not exclude on basis of low cost: Hence
                                % isnew1||isnew2||thisTensorXflag>0
                                [newCost isOK] = getBuildCost(freelegs,commonlegs,legCosts,costType,oldmuCap,muCap,isnew1||isnew2||thisTensorXflag>0,costToBuild1,costToBuild2);
                                if ~isOK
                                    newmuCap = min([newmuCap newCost]);
                                end
                            end
                            
                            % If cost is OK, compare with previous best known cost for constructing this object
                            if isOK
                                % Get involved tensors
                                tensorsInNew = tensorsIn1 | tensorsIn2;
                                % Find if previously constructed
                                objectPtr = 0;
                                for x=1:numel(objects{numInObjects})
                                    if isequal(objects{numInObjects}{x}{2},tensorsInNew)
                                        objectPtr = x;
                                        break;
                                    end
                                end
                                isnew = objectPtr==0;
                                if isnew
                                    % Is a new construction
                                    objectPtr = numel(objects{numInObjects})+1;
                                else
                                    % Compare new cost with best-so-far cost for construction of this object
                                    isOK = isLessThan(newCost,objects{numInObjects}{objectPtr}{4},costType);
                                end
                                
                                % ### If appropriate, update tensorXlist (list of tensors which can be contracted with objects created by outer product)
                                if allowOPs
                                    if isOK
                                        % New tensor or new best cost
                                        E_is_2 = ~any(freelegs2) && any(freelegs1);
                                        if (~any(freelegs1) && any(freelegs2)) || E_is_2
                                            % New best sequence consistent with Fig.5(c).
                                            % Determine the value of allIn, which corresponds to xi_C. (This is used in determining valid tensors X to contract with outer products).
                                            if E_is_2
                                                allIn = getProdLegDims(legs2,legCosts,costType);
                                            else
                                                allIn = getProdLegDims(legs1,legCosts,costType);
                                            end
                                            % Add to tensor X list for outer products (or if already there, update the value of allIn):
                                            if isnew
                                                if isGreaterThan_sd(dimSquared(allIn,costType),getProdLegDims(freelegs,legCosts,costType),costType) % Enforce Eq.(27)
                                                    [tensorXlegs tensorXdims tensorXflags] = addToTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs,allIn,costType,isnew);
                                                end
                                            else
                                                if isGreaterThan_sd(dimSquared(allIn,costType),getProdLegDims(freelegs,legCosts,costType),costType) % Enforce Eq.(27)
                                                    [tensorXlegs tensorXdims tensorXflags] = addToTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs,allIn,costType,isnew,~isequal(objects{numInObjects}{objectPtr}{7},zeros(1,costType)));
                                                elseif ~isequal(objects{numInObjects}{objectPtr}{7},zeros(1,costType)) % Only need to invoke removeFromTensorXList if there might actually be an entry
                                                    [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs);
                                                end
                                            end
                                        else
                                            % This tensor is not an eligible tensor X for an outer product: Store a dummy value in allIn to indicate this
                                            allIn = zeros(1,costType);
                                            % Best cost and not consistent with Fig.5(c): Ensure does not appear in tensorXlist. Active removal only
                                            % required if object is not new, and previous best sequence was consistent with Fig.5(c), so allIn is not a
                                            % dummy on the old entry.
                                            if ~isnew && ~isequal(objects{numInObjects}{objectPtr}{7},zeros(1,costType))
                                                [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs);
                                            end
                                        end
                                    elseif isequal(newCost,objects{numInObjects}{objectPtr}{4})
                                        % Equal-best cost to a known sequence for the same tensor
                                        if ~isequal(objects{numInObjects}{objectPtr}{7},zeros(1,costType))
                                            % Previous best sequence was consistent with Fig.5(c) so tensor may appear in the provisional environments list
                                            E_is_2 = ~any(freelegs2) && any(freelegs1);
                                            if (~any(freelegs1) && any(freelegs2)) || E_is_2
                                                % Determine the value of allIn, which corresponds to xi_C in Fig.5(c).
                                                if E_is_2
                                                    allIn = getProdLegDims(legs2,legCosts,costType);
                                                else
                                                    allIn = getProdLegDims(legs1,legCosts,costType);
                                                end
                                                % If smaller than previous value, update the value of allIn:
                                                if isGreaterThan_sd(objects{numInObjects}{objectPtr}{7},allIn,costType)
                                                    if isGreaterThan_sd(dimSquared(allIn,costType),getProdLegDims(freelegs,legCosts,costType),costType) % Enforce Eq.(27)
                                                        [tensorXlegs tensorXdims tensorXflags] = updateTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs,allIn);
                                                    else
                                                        [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs);
                                                    end
                                                end
                                            else
                                                % Found best-equal sequence not consistent with Fig.5(c)
                                                [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs);
                                            end
                                            %else: There already exists a best-known-cost sequence for the tensor which is not consistent with Fig.5(c), and isOK=false. Tensor does not appear in tensorXlist. No need to assign allIn. Now returning to start of main loop.
                                        end
                                        % else: Sequence is not capable of updating tensorXlist (not better cost, not equal cost). Also, isOK=false. No need to assign allIn. Now returning to start of main loop.
                                    end
                                else
                                    % Not doing outer products. Store a dummy value in allIn, which is never used.
                                    allIn = [];
                                end
                                % ### Done updating tensorXlist (list of tensors which can be contracted with objects created by outer product)
                            end
                            
                            if isOK
                                % Either no previous construction, or this one is better
                                if ~any(commonlegs)
                                    % ### This construction is an outer product. Note dimension of larger of the two participating tensors in newmaxdim. (This is used in enforcing index-dimension-related constraints.)
                                    newmaxdim = getProdLegDims(legs1,legCosts,costType);
                                    newmaxdim2 = getProdLegDims(legs2,legCosts,costType);
                                    if isGreaterThan_sd(newmaxdim2,newmaxdim,costType)
                                        newmaxdim = newmaxdim2;
                                    end
                                    % ### End recording dimension of larger of the two participating tensors in newmaxdim.
                                    thisIsOP = true;
                                    newseq = [seq1 seq2 0];
                                else
                                    % This construction is not an outer product.
                                    if isOP1
                                        newseq = [seq2 seq1 find(commonlegs)];
                                    else
                                        newseq = [seq1 seq2 find(commonlegs)];
                                    end
                                    thisIsOP = false;
                                    % ### This construction is not an outer product. Therefore store a dummy value in maxdim. (For outer products, maxdim records the dimension of the larger participating tensor, to assist in enforcing index-dimension-related constraints.)
                                    newmaxdim = [];
                                    % ### End storing dummy value in maxdim
                                end
                                
                                % Update objects{} with this construction
                                objects{numInObjects}{objectPtr} = {freelegs,tensorsInNew,newseq,newCost,thisIsOP,newmaxdim,allIn};
                                % ### Note 1: If this tensor has the structure of Fig.5(c) and so is capable of being contracted with an outer product object, |E| is recorded in allIn (otherwise this is a dummy value).
                                % ### Note 2: If this tensor is constructed by outer product, the dimension of the largest participating tensor is recorded in newmaxdim. (This is used in enforcing index-dimension-related constraints.) Otherwise, this is a dummy value.
                                
                                % Flag as a new construction
                                newobjectflags{numInObjects}(objectPtr) = true;
                                
                                % If top level, display result
                                if numInObjects == numtensors
                                    % ### If a valid contraction sequence has been found, there is no need to perform any contraction sequence more expensive than this. Set muCap accordingly.
                                    if costType==1
                                        muCap = newCost;
                                    else
                                        muCap = numel(newCost)-1;
                                    end
                                    % ### Done setting muCap accordingly.
                                    if verbosity > 1
                                        displayInterimCostAndSequence(newCost,newseq,costType,posindices,tracedindices);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % ### Finished searching if an object has been constructed which contains all tensors, and no new outer products have been enabled on the last pass (all(tensorXflags<2)==true, indicating that there are no new entries in the list of tensors which can be contracted with outer products).
    done = numel(objects{numtensors})~=0 && (all(tensorXflags<2) || ~allowOPs); % Final object has been constructed, and if outer products are allowed, also no new X's have recently been constructed
    
    if ~done
        if all(tensorXflags<2)
            % ### All X tensors have been present for an entire pass, so all permitted outer products at this cost have already been constructed.
            % Increment muCap, update oldmuCap
            if costType==1
                if newmuCap < muCap * max([min(legCosts) 2])
                    newmuCap = muCap * max([min(legCosts) 2]);
                end
            end
            oldmuCap = muCap;
            muCap = newmuCap;
            newmuCap = Inf;
        else
            % ### New X tensors generated this pass (some tensor X flags==2). Do another pass with same cost limit, to construct newly-allowed objects (i.e. sequences of affordable cost including at least one new outer product).
            % ### This is achieved by updating oldmuCap only. Now only outer products and contractions involving newly-created tensors will satisfy mu_0 < mu <= muCap.
            oldmuCap = muCap;
        end
        % Clear all new object flags
        for a=1:numel(newobjectflags)
            newobjectflags{a} = false(1,numel(newobjectflags{a}));
        end
        % ### Update tensor X flags (2 -> 1 -> 0):
        % ### 2: Newly created this pass becomes
        % ### 1: Created last pass; allow construction of cheap objects which contract with this, as they may have previously been excluded due to lack of a valid tensor X. 1 becomes...
        % ### 0: Old tensor X. Standard costing rules apply.
        % ### Delete redundant entries in tensorXlist (e.g. if A has a subset of the legs on B, and an equal or lower value of allIn (i.e. |E| in Fig.5(c)))
        [tensorXlegs tensorXdims tensorXflags] = mergeTensorXlist(tensorXlegs,tensorXdims,tensorXflags,costType);
        % ### Done updating tensor X flags
    else
        if numel(objects{numtensors})~=1
            error('Inappropriate number of objects containing all tensors!')
        end
    end
end

% Extract final result
sequence = int32(objects{numtensors}{1}{3});
cost = objects{numtensors}{1}{4};
end

function [newCost isOK] = getBuildCost(freelegs,commonlegs,legCosts,costType,oldmuCap,muCap,isnew,costToBuild1,costToBuild2)
% Get fusion cost
allLegs = freelegs | commonlegs;
isOK = true;
if costType==1
    % Is cost too high (>muCap)?
    newCost = prod(legCosts(allLegs)) + costToBuild1 + costToBuild2;
    if newCost > muCap
        isOK = false;
    end
    
    % Is cost too low (not made from new objects, and <=oldmuCap: This construction has been done before)
    if isOK && ~isnew
        if newCost <= oldmuCap
            isOK = false;
            newCost = Inf;
        end
    end
else
    % Is cost too high (>muCap)?
    fusionpower = sum(legCosts(allLegs,2));
    if fusionpower > muCap
        isOK = false;
        newCost = fusionpower;
    end
    
    % Is cost too low (not made from new objects, and <=oldmuCap: This construction has been done before)
    if isOK && ~isnew
        if fusionpower <= oldmuCap
            isOK = false;
            newCost = Inf;
        end
    end
    
    % If cost OK, determine total cost of construction
    if isOK
        newCostLen = max([numel(costToBuild1) numel(costToBuild2) fusionpower+1]);
        newCost = zeros(1,newCostLen);
        newCost(1:numel(costToBuild1)) = costToBuild1;
        if ~isempty(costToBuild2)
            newCost(1:numel(costToBuild2)) = newCost(1:numel(costToBuild2)) + costToBuild2;
        end
        newCost(fusionpower+1) = newCost(fusionpower+1) + prod(legCosts(allLegs,1));
    end
end
end

function flag = isLessThan(cost1,cost2,costType)
% Compares two full network costs
if costType==1
    flag = cost1 < cost2;
else
    if numel(cost1) < numel(cost2)
        flag = true;
    elseif numel(cost1) > numel(cost2)
        flag = false;
    else
        flag = false;
        for a = numel(cost2):-1:1
            if cost1(a) < cost2(a)
                flag = true;
                break;
            elseif cost1(a) > cost2(a)
                break;
            end
        end
    end
end
end

function flag = isGreaterThan_sd(cost1,cost2,costType)
% Compares two single-index costs
if costType==1
    flag = cost1 > cost2;
else
    flag = false;
    if cost1(2) > cost2(2)
        flag = true;
    elseif cost1(2) == cost2(2)
        if cost1(1) > cost2(1)
            flag = true;
        end
    end
end
end

function dim = getProdLegDims(freelegs,legCosts,costType)
if costType==1
    dim = prod(legCosts(freelegs));
else
    dim = [prod(legCosts(freelegs,1)) sum(legCosts(freelegs,2))];
end
end

function dim = dimSquared(dim,costType)
if costType==1
    dim = dim*dim;
else
    dim = [dim(1)*dim(1) dim(2)*2];
end
end

function displayInterimCostAndSequence(cost,seq,costType,posindices,tracedindices)
% Displays cost and sequence
dispseq = [tracedindices seq];
dispseq(dispseq>0) = posindices(dispseq(dispseq>0));
disp(['Sequence:       ' unpaddednum2str(dispseq)]);
if costType==1
    t = ['Cost:           ' num2str(cost)];
else
    t = 'Cost:           ';
    for a = numel(cost):-1:2
        t = [t num2str(cost(a)) 'X^' num2str(a-1) ' + ']; %#ok<AGROW>
    end
    t = [t num2str(cost(1)) 'X^0'];
end
if ~isempty(tracedindices)
    t = [t ' + tracing costs'];
end
disp(t);
end

function [tensorXlegs tensorXdims tensorXflags] = addToTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs,allIn,costType,isnew,oldMayHaveEntry)
% Constructed a new tensor for tensorXlist, or a known tensor at same or better cost.
% If legs exactly match a known non-provisional entry, consider as a possible tighter bound on allIn.
% Otherwise, consider for provisional list if not made redundant by any non-provisional entries.
% Add to provisional list if not made redundant by any non-provisional entries.
consider = true;
ptr = find(tensorXflags==1,1,'last');
for a=ptr:-1:1
    if isequal(tensorXlegs{a},freelegs)
        % If legs exactly match a non-provisional entry, update value of allIn.
        % If allIn for this entry just got increased, associated constraints have been relaxed. Flag this updated entry as provisional to trigger another pass.
        if isGreaterThan_sd(allIn,tensorXdims{a},costType)
            tensorXdims{a} = allIn;
            tensorXflags(a) = 2;
            % Flags are always in ascending order: Move re-flagged entry to end of list
            tensorXdims = tensorXdims([1:a-1 a+1:end a]);
            tensorXlegs = tensorXlegs([1:a-1 a+1:end a]);
            tensorXflags = tensorXflags([1:a-1 a+1:end a]);
        else
            tensorXdims{a} = allIn;
        end
        consider = false;
        break;
    else
        % Check to see if made redundant by existing non-provisional entry
        if all(tensorXlegs{a}(freelegs)) && ~isGreaterThan_sd(allIn,tensorXdims{a},costType)
            % All legs in freelegs are in overlap with tensorXlegs{a}, and dimension of absorbed tensor is not greater: Proposed new entry is redundant.
            % (Greater dimension is allowed as it would mean more permissive bounds for a subset of legs)
            consider = false;
            if ~isnew
                if oldMayHaveEntry
                    % This tensor: Excluded from list.
                    % Previous, higer-cost contraction sequence may have successfully made an entry. Remove it.
                    [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs);
                end
            end
            break;
        end
    end
end
if consider
    % Tensor not excluded after comparison with non-provisional entries:
    % If legs exactly match another provisional entry, new submission is a cheaper sequence so keep only the new value of allIn.
    % (This is the same subnetwork but with a better cost, and possibly a different value of \xi_C.)
    if ~isnew
        if oldMayHaveEntry
            for a=ptr+1:numel(tensorXlegs)
                if isequal(tensorXlegs{a},freelegs)
                    tensorXdims{a} = allIn;
                    consider = false;
                    break;
                end
            end
        end
    end
    % Otherwise: This is a new tensorXlist entry
    if consider
        tensorXlegs{end+1} = freelegs;
        tensorXflags(end+1) = 2;
        tensorXdims{end+1} = allIn;
    end
end
end

function [tensorXlegs tensorXdims tensorXflags] = updateTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs,allIn)
% Found new contraction sequence for an existing entry, but with a smaller value of allIn. Update the stored value to this value.
% (This is the same subnetwork, but the \xi_C constraint has been tightened.)
for a=numel(tensorXlegs):-1:1
    if isequal(freelegs,tensorXlegs{a})
        tensorXdims{a} = allIn;
        break;
    end
end
% If the original entry was provisional and redundant, this one is too. No match will be found (as the redundant tensor was never recorded), and that's OK.
% If the original entry was not provisional and redundant, it has now been updated.
end

function [tensorXlegs tensorXdims tensorXflags] = removeFromTensorXlist(tensorXlegs,tensorXdims,tensorXflags,freelegs)
% Remove this tensor from tensorXlist.
% Cannot restrict checking just to provisional portion, because might need to remove a confirmed tensor's data from the list in the scenario described in footnote 2.
for a=numel(tensorXlegs):-1:1
    if isequal(freelegs,tensorXlegs{a})
        tensorXlegs(a) = [];
        tensorXflags(a) = [];
        tensorXdims(a) = [];
        break;
    end
end
end

function [tensorXlegs tensorXdims tensorXflags] = mergeTensorXlist(tensorXlegs,tensorXdims,tensorXflags,costType)
% Merge provisional tensors into main list and decrease all nonzero flags by one.
% For each provisional entry in turn, check if it makes redundant or is made redundant by any other provisional entries, and if it makes redundant any
% non-provisional entries. If so, delete the redundant entries.

% Decrease non-zero tensorXflags
tensorXflags(tensorXflags>0) = tensorXflags(tensorXflags>0)-1;

ptr = find(tensorXflags==1,1,'first');
for a=numel(tensorXlegs):-1:ptr
    % Permit provisional tensors to be eliminated by other provisional tensors (elimination by non-provisional tensors was done earlier)
    for b=[ptr:a-1 a+1:numel(tensorXlegs)]
        if all(tensorXlegs{b}(tensorXlegs{a})) && ~isGreaterThan_sd(tensorXdims{a},tensorXdims{b},costType) % All legs in tensorXlegs{a} overlap with tensorXlegs{b}, and dimension of absorbed tensor is not greater: tensorXlegs{a} is redundant.
            tensorXlegs(a) = [];
            tensorXflags(a) = [];
            tensorXdims(a) = [];
            break;
        end
    end
end

for a=ptr-1:-1:1
    % Now permit non-provisional tensors to be eliminated by provisional tensors, as these are now confirmed
    for b=ptr:numel(tensorXlegs)
        if all(tensorXlegs{b}(tensorXlegs{a})) && ~isGreaterThan_sd(tensorXdims{a},tensorXdims{b},costType) % All legs in tensorXlegs{a} overlap with tensorXlegs{b}, and dimension of absorbed tensor is not greater: tensorXlegs{a} is redundant.
            tensorXlegs(a) = [];
            tensorXflags(a) = [];
            tensorXdims(a) = [];
            break;
        end
    end
end
end
