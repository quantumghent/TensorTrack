function wig = Wigner6j(j1,j2,j3,j4,j5,j6,ifs,ifcb)
% Computes the Wigner $6j$ symbols
%
% .. math::
%
%    W = \begin{Bmatrix}
%       j_1 & j_2 & j_3 \\
%       j_4 & j_5 & j_6
%    \end{Bmatrix}
% 
% or the recoupling matrix element
%
% .. math::
%
%    W = \left\langle (j_1,j_2)j_3,j_4,j_5 \middle| j_1,(j_2,j_4)j_6,j_5 \right\rangle
%
% using the Racah formula.
%
% VERSION of January, 2020
%
% Programmer: V. B. Sovkov, St. Petersburg State University, Shanxi University
%
% Usage
% -----
% :code:`W = Wigner6j(j1,j2,j3,j4,j5,j6,ifs,ifcb)`
% :code:`Wigner6j(j1,j2,j3,j4,j5,j6)`
%
% Arguments
% ---------
% j1, j2, j3, j4, j5, j6
%   set of six angular momenta, must be arrays of the same size
% 
% Optional Arguments
% ------------------
% ifs
%   the switch of the computational methods;
%   although the central part of all the computations is based on the Racah formula, some details differ:
%
% 	- 0 - the symbolic (presumably accurate) computation with the double-precision output
%  	- -1 - the same symbolic computation as ifs=0 but with the symbolic output (accurate
%     square root/rational-type equation, simplified)
% 	- -2 - the same as :code:`ifs=-1` but without a final simplification of the symbolic
%     expression, which can be time-consuming sometimes
%  	- 1, 2 (default, recommended), 3 - numeric double-precision algorithms (see Notes below)
%
%     all other input values of ifs are set to the closest from the above list.
% 
% ifcb
% 	if exists and is true, switches to computing the coupling matrix elements instead of the
%   $6j$-symbols (default :code:`false`)
%
% Returns
% -------
% W
%   the resulting values of the 3j-symbols (:code:`ifcb=false`) or the Clebsch-Gordan
%   coefficients (:code:`ifcb=true`) in either numeric double-precision or symbolic form
%   (see the input parameter :code:`ifs`); array of the same size as $j_1$.
% 
% Notes
% -----
% most of the estimates below is based on extensive numerical tests within a range of the
% quantum nubers <=1000
%
% a. the numeric algorithms (:code:`ifs>0`) are usually much faster than the symbolic
%    algorithm (:code:`ifs<=0`)
% b. the accuracy of the symbolic algorithm remains the uttermost possible
% c. the accuracy of the numeric algorithms can worsen for big quantum numbers
% d. all the "numeric" algorithms switch automatically to the symbolic computations
% e. as soon as the numeric overflow occurs, thereby improving the accuracy
%    but slowing down the computations for some big quantum numbers
% f. in cases with at least one of the $j$ quantum numbers <=2, the explicit accurate
%    equations are applied in all the algorithms ensuring the highest possible accuracy
% g. for relatively small quantum numbers (up to ~20) all the algorithms provide
%    approximately equivalent results with a reasonably high accuracy
% h. the algorithms :code:`ifs=1` and :code:`ifs=2` ensure a reasonably good accuracy; the
%    worst registered ca for both algorithms (occuring before switching to the symbolic
%    regime) was :code:`Wigner6j(41.5,52,43.5,36,38.5,46)=-0.00029` with the biggest
%    relative inaccuracy of 3.6e-10
% i. the algorithm :code:`ifs=2` (comparing to :code:`ifs=1`) proceeds in a wider range of
%    quantum numbers numerically before switching to the symbolic regime (i.e., can work
%    faster with big quantum numbers); in our experimental probing it proved to ensure a
%    good enough accuracy and was chosen as a default one
% j. the algorithm :code:`ifs=3` (comparing to :code:`ifs=2`) proceeds in a wider range of
%    quantum numbers numerically before switching to the symbolic regime (i.e., can work
%    faster with big quantum numbers) but for some combinations of big quantum numbers,
%    completely wrong results are observed: e.g.,
%    :code:`Wigner6j(227,230.5,210.5,249.5,235,232,3)`; thereby, we do not recommend using
%    the algorithm code:`ifs=3` at all, and only keep it here in view of the completeness.
%
% References
% ----------
% 1. https://en.wikipedia.org/wiki/6-j_symbol
% 2. `Zare, R. N., & Harter, W. G. (1988). Angular momentum: understanding spatial aspects
%    in chemistry and physics. New York, 120.
%    <http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471858927.html>`_

NJ = max([numel(j1), numel(j2), numel(j3), numel(j4), numel(j5), numel(j6)]);

if numel(j1) ~= NJ && isscalar(j1)
    j1(1:NJ) = j1;
end
if numel(j2) ~= NJ && isscalar(j2)
    j2(1:NJ) = j2;
end
if numel(j3) ~= NJ && isscalar(j3)
    j3(1:NJ) = j3;
end
if numel(j4) ~= NJ && isscalar(j4)
    j4(1:NJ) = j4;
end
if numel(j5) ~= NJ && isscalar(j5)
    j5(1:NJ) = j5;
end
if numel(j6) ~= NJ && isscalar(j6)
    j6(1:NJ) = j6;
end

try


    NJ = min([numel(j1) numel(j2) numel(j3) numel(j4) numel(j5) numel(j6)]);
    Nj1S=size(j1);
    j1 = double(reshape(j1(1:NJ),NJ,1));
    j2 = double(reshape(j2(1:NJ),NJ,1));
    j3 = double(reshape(j3(1:NJ),NJ,1));
    j4 = double(reshape(j4(1:NJ),NJ,1));
    j5 = double(reshape(j5(1:NJ),NJ,1));
    j6 = double(reshape(j6(1:NJ),NJ,1));

    if nargin<7 || isempty(ifs) || ~isnumeric(ifs) && ~islogical(ifs) || ifs(1)>0 && ifs(1)<=1.5 % choose the algorithm
        ifs = 2; % default
    elseif ifs(1)>0
        ifs = min(round(ifs(1)),3);
    else
        ifs=max(round(ifs(1)),-2);
    end

    NK = if3jc(j1,j2,j3) & if3jc(j1,j5,j6) & if3jc(j3,j4,j5) & if3jc(j2,j4,j6);
    if all(~NK)
        if ifs>=0
            wig=zeros(NJ,1);
        else
            wig=sym(zeros(NJ,1));
        end
        if prod(Nj1S)==NJ
            reshape(wig, Nj1S);
        end
        return;
    end

    if ifs>0
        wig=zeros(NJ,1);
    else
        wig=sym(zeros(NJ,1));
    end

    ifcb = nargin>7 && ~isempty(ifcb) && ifcb(1); % if the Clebsch-Gordan coefficient instead of the 3j-symbols are needed?
    if ifcb
        J1=j1;
        J2=j2;
        J3=j3;
        J4=j4;
        J5=j5;
        J6=j6;
    end

    for k=1:NJ
        jup=[j1(k),j2(k),j3(k)]';
        jdw=[j4(k),j5(k),j6(k)]';
        if min(jdw)<min(jup)
            [~,jup,jdw]=ArranA(jdw,jup);
            jj=jdw(3);
            jdw(3)=jup(3);
            jup(3)=jj;
        end
        [~,jup,jdw]=ArranA(jup,jdw);
        while jup(2)>jdw(2)
            jj=jdw(2:3);
            jdw(2:3)=jup(2:3);
            jup(2:3)=jj;
            [~,jup,jdw]=ArranA(jup,jdw);
        end
        j1(k)=jup(1);
        j2(k)=jup(2);
        j3(k)=jup(3);
        j4(k)=jdw(1);
        j5(k)=jdw(2);
        j6(k)=jdw(3);
    end
    if NJ>1
        [ind0,j1,j2,j3,j4,j5,j6] = ...
            ArranA ( reshape(j1(1:NJ),NJ,1),reshape(j2(1:NJ),NJ,1),reshape(j3(1:NJ),NJ,1),reshape(j4(1:NJ),NJ,1),reshape(j5(1:NJ),NJ,1),reshape(j6(1:NJ),NJ,1));
        NK = NK(ind0);

        kk = find ( j1(1:NJ-1)==j1(2:NJ) & j2(1:NJ-1)==j2(2:NJ) & j3(1:NJ-1)==j3(2:NJ) & j4(1:NJ-1)==j4(2:NJ) & j5(1:NJ-1)==j5(2:NJ) & j6(1:NJ-1)==j6(2:NJ) );
        if ~isempty(kk)
            NK(kk)=0;
        end
    else
        kk=[];
    end

    NKs=[]; % for pointing to erronious numeric results in order to switch them to the symbolic computations

%%
%   special cases - faster and more accurate computation than the general case below

    if any(NK & j1<=2) % the range of the special cases currently included

%   case j = 0
        kp = find(NK & j1==0);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = (-1).^(j2(kp)+j4(kp)+j5(kp)) ./ sqrt((2*j2(kp)+1).*(2*j5(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^(j2(kp)+j4(kp)+j5(kp)) ./ sqrt(sym((2*j2(kp)+1).*(2*j5(kp)+1)));
            end
            NK(kp)=0;
        end
%   case   j = 1/2
        kp = find(NK & j1==1/2 & j5<j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s-2*j5(kp))./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j3(kp)+1)./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s-2*j5(kp))./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j3(kp)+1)./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1/2); %  & j5>j6
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
%   case   j = 1
        kp = find(NK & j1==1 & j5>j6 & j3>j2); % & j3>j2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(s./(2*j5(kp)-1)./(2*j5(kp)).*(s+1)./(2*j5(kp)+1).*(s-2*j4(kp)-1)./(2*j3(kp)-1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(s./(2*j5(kp)-1)./(2*j5(kp)).*(s+1)./(2*j5(kp)+1).*(s-2*j4(kp)-1)./(2*j3(kp)-1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5<j6 & j3>j2); % & j3>j2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s-2*j5(kp)-1)./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j5(kp))./(2*j5(kp)+3)./(2*j3(kp)-1).*(s-2*j3(kp)+1)./(2*j3(kp)).*(s-2*j3(kp)+2)./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s-2*j5(kp)-1)./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j5(kp))./(2*j5(kp)+3)./(2*j3(kp)-1).*(s-2*j3(kp)+1)./(2*j3(kp)).*(s-2*j3(kp)+2)./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5==j6 & j3>j2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(2*(s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j5(kp)+2)./(2*j3(kp)-1).*(s-2*j5(kp))./(2*j3(kp)).*(s-2*j3(kp)+1)./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(2*(s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j5(kp)+2)./(2*j3(kp)-1).*(s-2*j5(kp))./(2*j3(kp)).*(s-2*j3(kp)+1)./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5<j6 & j3==j2);
        if ~isempty(kp)
            s=j2(kp)+j4(kp)+j6(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(2*(s+1)./(2*j2(kp))./(2*j2(kp)+1).*(s-2*j4(kp))./(2*j2(kp)+2)./(2*j6(kp)-1).*(s-2*j2(kp))./(2*j6(kp)).*(s-2*j6(kp)+1)./(2*j6(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(2*(s+1)./(2*j2(kp))./(2*j2(kp)+1).*(s-2*j4(kp))./(2*j2(kp)+2)./(2*j6(kp)-1).*(s-2*j2(kp))./(2*j6(kp)).*(s-2*j6(kp)+1)./(2*j6(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5==j6  & j2==j3); % & j2==j3
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = 2*(-1).^s .* (j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1)) ./ sqrt((2*j5(kp)).*(2*j5(kp)+1).*(2*j5(kp)+2).*(2*j3(kp)).*(2*j3(kp)+1).*(2*j3(kp)+2));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .* sym(j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1)) ./ sqrt(sym((2*j5(kp)).*(2*j5(kp)+1).*(2*j5(kp)+2).*(2*j3(kp)).*(2*j3(kp)+1).*(2*j3(kp)+2)));
            end
            NK(kp)=0;
        end
%   case   j = 3/2
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==3/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+1;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       (s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       (s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==1/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==3/2);
        if ~isempty(kp)
            s=j4(kp)+j3(kp)+j5(kp);
			A1=2*j3(kp)+2;
			A2=2*j5(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==-3/2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==-1/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==-3/2); % & j5-j6==-3/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       (s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       (s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j5-j6==1/2 & j3-j2==1/2); % & j3-j2==1/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    ) .* ...
                       (2*(s-2*j5(kp)).*(s-2*j3(kp))-(s+2).*(s-2*j4(kp)-1))...
				;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    )) .* ...
                       sym(2*(s-2*j5(kp)).*(s-2*j3(kp))-(s+2).*(s-2*j4(kp)-1))...
				;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==-1/2); % & j3-j2==1/2 & j5-j6==-1/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    ) .* ...
                       ((s-2*j5(kp)-1).*(s-2*j3(kp))-2*(s+2).*(s-2*j4(kp)))...
				;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    )) .* ...
                       sym((s-2*j5(kp)-1).*(s-2*j3(kp))-2*(s+2).*(s-2*j4(kp)))...
				;
            end
            NK(kp)=0;
        end
%   case   j = 2, j3-j2=2
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+1;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2).*(s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-3).*(s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2).*(s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-3).*(s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5==j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==-1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==-2); % & j5-j6==-2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+5;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)-3).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3).*(s-2*j3(kp)+4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)-3).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3).*(s-2*j3(kp)+4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
%   case   j = 2, j3-j2=1
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==2);
        if ~isempty(kp)
            s=j4(kp)+j3(kp)+j5(kp);
			A1=2*j3(kp)+2;
			A2=2*j5(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==-2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+4;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-2).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2).*(s-2*j6(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-2).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2).*(s-2*j6(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 4*(-1).^s .*...
                    sqrt(...
                      s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)).*(j4(kp)-j5(kp)+1)-(j3(kp)-1).*(j3(kp)-j5(kp)+1) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 4*(-1).^s .*...
                    sqrt(sym(...
                      s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)).*(j4(kp)-j5(kp)+1)-(j3(kp)-1).*(j3(kp)-j5(kp)+1) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5==j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)+1).*(j4(kp)-j5(kp))-j3(kp).^2+1 )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)+1).*(j4(kp)-j5(kp))-j3(kp).^2+1 )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==-1); % & j5-j6==-1
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 4*(-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)+2).*(j4(kp)-j5(kp)-1)-(j3(kp)-1).*(j5(kp)+j3(kp)+2) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 4*(-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)+2).*(j4(kp)-j5(kp)-1)-(j3(kp)-1).*(j5(kp)+j3(kp)+2) )...
                    ;
            end
            NK(kp)=0;
        end
%   case   j = 2, j3==j2
        kp = find(NK & j1==2 & j3==j2 & j5-j6==-2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)-1).*(s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)-1).*(s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3==j2 & j5-j6==-1);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+2;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)).*(s-2*j6(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j2(kp)+1).*(j4(kp)-j2(kp))-j6(kp).^2+1 )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)).*(s-2*j6(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j2(kp)+1).*(j4(kp)-j2(kp))-j6(kp).^2+1 )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3==j2 & j5==j6); % & j3==j2 & j5==j6
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+3;
			C =j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1);
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    (3*C.*(C+1)-4*j5(kp).*(j5(kp)+1).*j3(kp).*(j3(kp)+1)) ./...
                    sqrt(...
                    A1.*(A1-1).*(A1-2).*(A1-3).*(A1-4)...
                    .* A2.*(A2-1).*(A2-2).*(A2-3).*(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sym(3*C.*(C+1)-4*j5(kp).*(j5(kp)+1).*j3(kp).*(j3(kp)+1)) ./...
                    sqrt(sym(...
                    A1.*(A1-1).*(A1-2).*(A1-3).*(A1-4)...
                    .* A2.*(A2-1).*(A2-2).*(A2-3).*(A2-4)...
				    ));
            end
            NK(kp)=0;
        end

    end
%%

    NK = find(NK);

    if ~isempty(NK)
        if ifs==3
            C =     tc(j1(NK),j2(NK),j3(NK),ifs)  + tc(j1(NK),j5(NK),j6(NK),ifs)  + tc(j4(NK),j2(NK),j6(NK),ifs)  + tc(j4(NK),j5(NK),j3(NK),ifs);
        else
            C =     tc(j1(NK),j2(NK),j3(NK),ifs) .* tc(j1(NK),j5(NK),j6(NK),ifs) .* tc(j4(NK),j2(NK),j6(NK),ifs) .* tc(j4(NK),j5(NK),j3(NK),ifs);
        end
        k0 = find(isinf(C) | isnan(C));
        if ~isempty(k0)
            NKs = [NKs;NK(k0)];
            NK(k0) = [];
            C(k0)=[];
        end
        
        if ~isempty(NK)

            t1 = max ( [j1(NK) + j2(NK) + j3(NK), j1(NK) + j5(NK) + j6(NK), j4(NK) + j2(NK) + j6(NK), j4(NK) + j5(NK) + j3(NK)] , [] , 2 );
            t2 = min ( [j1(NK) + j2(NK) + j4(NK) + j5(NK), j2(NK) + j3(NK) + j5(NK) + j6(NK), j3(NK) + j1(NK) + j6(NK) + j4(NK)] , [] , 2 );

            if ifs>0 % numeric
                if ifs==3
                    for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                        kN = NK(k);
                        t = (t1(k):t2(k))';
                        wig(kN) = ((-1).^t)' * exp ( gammaln([...
                            (t+2),(t-j1(kN)-j2(kN)-j3(kN)+1),(t-j1(kN)-j5(kN)-j6(kN)+1),(t-j4(kN)-j2(kN)-j6(kN)+1),(t-j4(kN)-j5(kN)-j3(kN)+1),(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t+1),(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t+1),(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t+1)...
                            ]) * [1;-1;-1;-1;-1;-1;-1;-1] + C(k) );
                    end
                else % ifs==1 || ifs==2
                    for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                        kN = NK(k);
                        t = (t1(k):t2(k))';
                        wig(kN) = (-1).^t' * (...
                            factorial(t+1)...
                            ./factorial(t-j1(kN)-j2(kN)-j3(kN))...
                            * C(k)...
                            ./factorial(t-j1(kN)-j5(kN)-j6(kN))...
                            ./factorial(t-j4(kN)-j2(kN)-j6(kN))...
                            ./factorial(t-j4(kN)-j5(kN)-j3(kN))...
                            ./factorial(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t)...
                            ./factorial(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t)...
                            ./factorial(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t)...
                            );
                    end
                end
                k0 = find(isinf(wig(NK)) | isnan(wig(NK)));
                if ~isempty(k0)
                    NKs = [NKs;NK(k0)];
                end
            else % symbolic
                for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                    kN = NK(k);
                    t = (t1(k):t2(k))';
                    wig(kN) = (-1).^t' * ( ...
                        factorial(sym(t+1))...
                        ./factorial(sym(t-j1(kN)-j2(kN)-j3(kN)))...
                        * C(k)...
                        ./factorial(sym(t-j1(kN)-j5(kN)-j6(kN)))...
                        ./factorial(sym(t-j4(kN)-j2(kN)-j6(kN)))...
                        ./factorial(sym(t-j4(kN)-j5(kN)-j3(kN)))...
                        ./factorial(sym(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t))...
                        ./factorial(sym(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t))...
                        ./factorial(sym(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t))...
                        );
                end
            end
        end
    end
%%
%   post-processing

    if~ifcb
        if ~ifs
            wig = double(wig);
        elseif ifs==-1
            wig = simplify(wig);
        end
    end

    if ~isempty(NKs) % incorrect numerical results - swich to the symbolic computations
        NKs=sort(NKs);
        wig(NKs) = Wigner6j(j1(NKs),j2(NKs),j3(NKs),j4(NKs),j5(NKs),j6(NKs),0);
    end

    if ~isempty(kk)
        k0 = diff(kk);
        if all(k0~=1)
            wig(kk) = wig(kk+1);
        else
            k0 = find(k0>1);
            k00=1;
            for k=1:numel(k0)
                wig(kk(k00:k0(k))) = wig(kk(k0(k))+1);
                k00 = k0(k)+1;
            end
            wig(kk(k00:end)) = wig(kk(end)+1);
        end
     end

    if exist('ind0','var')
        wig(ind0)=wig;
    end

    if ifcb % recompute to the coupling matrix element
        if ifs>0
            wig = wig .* sqrt((2*J3+1).*(2*J6+1)) .* (-1).^(J1+J2+J4+J5);
        elseif ifs==0
            wig = double( wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5) );
        elseif ifs==-1
            wig = simplify( wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5) );
        else
            wig = wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5);
        end
    end
    
    if NJ>1 && prod(Nj1S)==NJ
        wig=reshape(wig,Nj1S);
    end

catch mectc
%%
    beep;
    disp(mectc);
    disp('Wigner6j ERROR: The input arguments can be incorrect; see comments in the code,');
    if ifs<=0 || exist('NKs','var') && ~isempty(NKs)
        disp('or the ``Symbolic'' toolbox is not properly installed')
    end
    try
        if NJ>1 && prod(Nj1S)==NJ
            wig=NaN(Nj1S);
        else
            wig=NaN(size(j1));
        end
    catch
        wig=NaN;
    end
    beep;
end

return;
end



function [ind0,varargout] = ArranA (varargin)
% simultaneous sorting of several input column vectors
% so that varargout{1} becomes an acsending-order version of varargin{1},
% and varargout{2} is resorted so that elements of varargin{2} corresponding
% to equal elements of i1 becomes an ascending-order subvectors as well, etc.;
% ind0  is a permutation vector.
%
% USAGE: [ind0,i1,i2,...] = ArranA (j1,j2,...)
% with j1, j2, ... being numerical arrays (matrices) containing the vectors to be arranged column-wise.
%
% Programmer: V. B. Sovkov
% St. Petersburg State University
% Shanxi University
%%
try
    if ~nargin
        ind0=[];
        varargout={};
        return;
    end
    [varargout{nargin},ind0] = sort(varargin{nargin},1);
    for m=1:size(ind0,2)
        for k=nargin-1:-1:1
            varargout{k}(:,m) = varargin{k}(ind0(:,m),m);
        end
    end
    for n=nargin-1:-1:1
        [varargout{n},ind] = sort(varargout{n},1);
        for m=1:size(ind0,2)
            for k=nargin:-1:1
                if k~=n
                    varargout{k}(:,m) = varargout{k}(ind(:,m),m);
                end
            end
            ind0(:,m)=ind0(ind(:,m),m);
        end
    end
catch mectc
    varargout=[];
    disp(mectc);
end



return
end



function if3j = if3jc(J1,J2,J3)
% checks necessary conditions on the angular momenta J1, J2, J3, which must be:
% (1) non-negative;
% (2) either integer or half-integer;
% (3) their sum must be integer;
% (4) triangular rule.
% returns true if all the conditions are fulfilled; otherwise returns false.

if3j = J1(:)>=0 & J2(:)>=0 & J3(:)>=0 & mod(J1(:)-J2(:)-J3(:),1)==0 & mod(2*J1(:),1)==0 & mod(2*J2(:),1)==0 & mod(2*J3(:),1)==0 & J3(:) <= J1(:)+J2(:) & J3(:) >= abs(J1(:)-J2(:));

return;
end







function tri = tc(j1,j2,j3,ifs)
% sqrt of the triangle coefficient or its logarithm (when ifs=3)
if nargin<4 || isempty(ifs) || ~isnumeric(ifs) || ~isreal(ifs) || ifs>0
    if ifs==3
        tri = zeros(size(j1));
        for k=1:numel(tri)
            tri(k) = gammaln(j1(k)+j2(k)-j3(k)+1)/2 + gammaln(j1(k)-j2(k)+j3(k)+1)/2 - sum( log( (-j1(k)+j2(k)+j3(k)+1):(j1(k)+j2(k)+j3(k)+1) )/2 );
        end
    else
        if ifs==1
            tri = sqrt( factorial(j1+j2+j3+1) ./ factorial(j1+j2-j3) );
        elseif ifs==2
            for k=numel(j1):-1:1
                tri(k,1) = prod( sqrt( (j1(k)+j2(k)-j3(k)+1):(j1(k)+j2(k)+j3(k)+1) ) );
            end
        end
        k=find(~isinf(tri) & ~isnan(tri));
        if ~isempty(k)
            tri(k) = sqrt( factorial(-j1(k)+j2(k)+j3(k)) ) ./ tri(k) .* sqrt( factorial(j1(k)-j2(k)+j3(k)) );
        end
    end
else
    tri = sqrt( factorial(sym(j1+j2-j3))./ factorial(sym(j1+j2+j3+1)) .* factorial(sym(j1-j2+j3)) .* factorial(sym(-j1+j2+j3)) );
end
return;
end
