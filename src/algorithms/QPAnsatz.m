classdef QPAnsatz
    % Quasi-Particle excitation ansatz
    
    properties
        alg_eigs = Arnoldi('MaxIter', 100, 'KrylovDim', 30, 'Tol', 1e-8, 'verbosity', Verbosity.diagnostics)
        alg_environments = struct('Tol', 1e-10, 'Algorithm', 'bicgstabl');
        howmany = 1
        which = 'smallestreal'
    end
    
    methods
        function alg = QPAnsatz(kwargs)
            arguments
                kwargs.?QPAnsatz
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    if isstruct(kwargs.(field{1}))
                        for field2 = fieldnames(kwargs.(field{1})).'
                            alg.(field{1}).(field2{1}) = kwargs.(field{1}).(field2{1});
                        end
                    else
                        alg.(field{1}) = kwargs.(field{1});
                    end
                end
            end
        end
        
        function [qp, mu] = excitations(alg, mpo, qp)
            if period(mpo) ~= period(qp)
                error('QPAnsatz:argerror', ...
                    'periodicity of mpo (%d) should be equal to that of the mps (%d)', ...
                    period(mpo), period(qp));
            end
            
            
            % Renormalization
            [GL, GR, lambda] = environments(mpo, qp.mpsleft, qp.mpsleft);
            
            for i = period(mpo):-1:1
                T = AC_hamiltonian(mpo, qp.mpsleft, GL, GR, i);
                offsetL(i) = dot(qp.mpsleft.AC{i}, apply(T{1}, qp.mpsleft.AC{i}));
            end
            [GL, GR, lambda] = environments(mpo, qp.mpsright, qp.mpsright);
            
            for i = period(mpo):-1:1
                T = AC_hamiltonian(mpo, qp.mpsright, GL, GR, i);
                offsetR(i) = dot(qp.mpsright.AC{i}, apply(T{1}, qp.mpsright.AC{i}));
            end
            offset = (offsetL + offsetR) / 2;
            
            GL = leftenvironment(mpo, qp.mpsleft, qp.mpsleft);
            
            % Algorithm
            H_effective = @(qp) updateX(alg, mpo, qp, GL, GR, offset);
            [qp, mu] = eigsolve(alg.alg_eigs, H_effective, qp, alg.howmany, alg.which);
            
            for i = alg.howmany:-1:2
                qp(i).X = X(i);
                qp(i).B = computeB(qp(i));
            end
        end
        
        function qp = updateX(alg, mpo, qp, GL, GR, offset)
            qp.B = computeB(qp);
            B = qp.B;
            
            H_c = B_hamiltonian(mpo, qp, GL, GR, 'Type', 'center');
            for i = period(qp):-1:1
                B{i} = MpsTensor(apply(H_c{i}, B{i}), 1);
            end
            
            H_l = B_hamiltonian(mpo, qp, GL, GR, 'Type', 'left');
            for i = 1:period(qp)
                B{i} = B{i} + repartition(apply(H_l{i}, qp.AR{i}), rank(B{i}));
            end
           
            H_r = B_hamiltonian(mpo, qp, GL, GR, 'Type', 'right');
            for i = 1:period(qp)
                B{i} = B{i} + repartition(apply(H_r{i}, qp.AL{i}), rank(B{i}));
            end
            
            for i = 1:period(qp)
                qp.B{i} = B{i} - qp.B{i} * offset(i);
            end

            qp.X = computeX(qp);
        end
    end
end

