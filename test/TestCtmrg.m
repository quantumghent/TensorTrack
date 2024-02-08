classdef TestCtmrg < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        pspace = struct('fermion', GradedSpace.new(fZ2(0, 1), [1 1], false))
        vspace = struct('fermion', GradedSpace.new(fZ2(0, 1), [1 1], false))
        chispace =  struct('fermion', GradedSpace.new(fZ2(0, 1), [10 10], false))
        alg = {Ctmrg('tol', 1e-6)}
        %other symms TBA
        %larger unit cells TBA
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, pspace, vspace, chispace)
            peps = UniformPeps(PepsTensor.randnc(pspace, vspace, vspace));
            envs = CtmrgEnvironment(peps, peps, chispace);
            %check matching spaces TBA
        end

        function testConvergence(tc, pspace, vspace, chispace, alg)
            maxloops = 10;
            i = 0;
            err = 1;
            while err > alg.tol && i < maxloops
                i = i+1;
                peps = UniformPeps(PepsTensor.randnc(pspace, vspace, vspace));
                envs = CtmrgEnvironment(peps, peps, chispace);
                [envs, norm, err] = fixedpoint(alg, peps, peps, envs);                
            end
            assert(i < maxloops, "Convergence seems to fail for simple random PEPS.")
            %check left move TBA
            %check norm TBA
        end

        function testPWave(tc)
            %initialize
            pspace = GradedSpace.new(fZ2(0, 1), [1 1], false);
            vspace = pspace;
            chispace =  GradedSpace.new(fZ2(0, 1), [10 10], false);
            alg = Ctmrg('tol', 1e-6, 'miniter', 1, 'maxiter', 200);
            %fill PEPS as PWave from Gaussian
            mblocks{1} = [-0.4084 - 0.5139i  -0.1033 + 0.1492i  -0.7413 - 1.1574i   0.5924 + 0.6552i ...
                          -0.0660 - 0.1002i  -0.1678 - 0.3237i  -0.3739 + 0.9139i   0.6497 + 0.8648i];
            mblocks{2} = [ 0.5454 + 0.6080i  -0.3003 - 0.2482i   0.3430 - 0.1734i  -0.2851 - 0.2658i ...
                          -0.4272 + 0.5267i   0.4567 + 0.5106i  -0.6708 - 0.1306i   0.2184 - 0.1094i];
            A = PepsTensor.zeros(pspace, vspace, vspace);
            peps = UniformPeps(PepsTensor(A.var.fill_matrix(mblocks)));
            A = peps.A{1}.var;
            %run ctmrg
            maxloops = 10;
            i = 0;
            err = 1;
            while err > alg.tol && i < maxloops
                i = i+1;
                envs = CtmrgEnvironment(peps, peps, chispace);
                [envs, norm, err] = fixedpoint(alg, peps, peps, envs);                
            end
            assert(i < maxloops, "Convergence seems to fail for PWave example.")
            %local density comparison to Gaussian
            O = Tensor(pspace, pspace);
            O.var.var{2} = 1;
            GL = contract(envs.corners{1,1,1}, [1 -4], envs.edges{4,1,1}, [2 -2 -3 1], envs.corners{4,1,1}, [-1 2]);
            GR = contract(envs.corners{2,1,1}, [-1 1], envs.edges{2,1,1}, [1 -2 -3 2], envs.corners{3,1,1}, [2 -4]);
            n = contract(GL, 1:4, GR, [4,2,3,1]);
            GL = GL/sqrt(n);
            GR = GR/sqrt(n);
            l = contract(envs.edges{3,1,1}, [-1 5 2 1], GL, [1 4 3 8], envs.edges{1,1,1}, [8 9 10 -4], ...
                    A, [7 4 5 -2 9], conj(A), [6 3 2 -3 10], O, [6 7]);
            lw = contract(envs.edges{3,1,1}, [-1 5 2 1], GL, [1 4 3 8], envs.edges{1,1,1}, [8 9 10 -4], ...
                    A, [6 4 5 -2 9], conj(A), [6 3 2 -3 10]);
            rho = contract(l, [2 3 4 1], GR, [1 3 4 2])/contract(lw, [2 3 4 1], GR, [1 3 4 2]);
            assert(abs(rho - 0.178787157813086) < sqrt(alg.tol), "CTMRG yields wrong local density for PWave example.")
        end
    end
end

