function funs = auto_acpc_reorient_autogen
% funs = auto_acpc_reorient_autogen
% Genetic algorithm based optimization functions for neuroimaging registration.
%
% License: MIT License
% Copyright (C) 2020-2024 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    funs = importFunctions; % For Octave compatibility: we need the first function to have the same name as the filename.
end

% For MatLab compatibility, we use a function to return other functions handlers as properties, this is a workaround since it cannot load multiple functions in one .m file (contrarywise to Octave using source())
% This way, we can just do the following to call any function here from any other script:
% to load the aux lib: addpath(genpath(strcat(cd(fileparts(mfilename('fullpath'))),'/../gbnn-core/'))); aux = gbnn_aux;
% to use a function: aux.func_name(args);
function funs = importFunctions
    funs.mga_imatrix = @mga_imatrix;
end  % endfunction

function [best, bestfit, bestfitstd, bestpool, bestpoolfit] = mga_imatrix(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, maxdim, randomgen, mutatefunc, truerandom, onlyrigid, debugmode)
    % [best, bestfit, bestfitstd] = mga_imatrix(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, maxdim, randomgen, mutatefunc, onlyrigid, debugmode)
    % Microbial Genetic Algorithm to optimize an spm_imatrix() vector. This will find a better/locally optimal translation and orientation given the provided fitnessfunc by using the microbial genetic algorithm (MGA). The advantage compared to using stochastic descent as SPM does is that genetic algorithms can explore the whole solution space, even if it is non-continuous or non-differentiable, and is not affected by the curse of dimensionality.
    %
    % Inputs:
    % * genestart: a vector of a single candidate that can be fed to fitnessfunc(genestart). This will both be used as the reference candidate (so we ensure we don't return a worse solution), and to generate new candidates. genestart should be a 12 elements vector as made by spm_imatrix(svol.mat), but don't forget that spm_coreg.optfun and spm_hist2 will apply it on top of input_vol.mat, such that orientation_that_will_be_tested = pinv(spm_matrix(genestart))*input_vol.mat -- in spm_coreg you can see the full call to spm_hist2 being: VF.mat\spm_matrix(x(:)')*VG.mat -- so by default you want genestart to be [0 0 0 0 0 0] (as is done in spm_coreg, which means we start from input_vol.mat as-is without any modification), NOT spm_imatrix(input_vol.mat) which would be a transform that would be applied ON TOP of input_vol.mat
    % * maxdim: maximum dimension to generate rotations when combined with randomstart = true. In general, set maxdim = svol.dim
    % * popsize: size of the population. This number should be relative to the number of tournamentcount, a good number is 1/10th of tournamentcount, as to allow all candidates to be explored and recombined by the end of the algorithm.
    % * demesize: size of the deme, which is the subpopulation size inside the ring that will be used to find a candidate B to compare candidate A. See Inman Harvey's paper. A good value is demesize = 0.1 to 0.3 * popsize.
    % * fitnessfunc: fitness function to use, it should accept a single variable x which will be a vector of the same size as genestart. Tip: try to optimize this function to be super fast, as although this genetic algorithm uses caching to try to save some time, a computation expensive fitnessfunc will limit the number of tournamentcount you can do (and hence the likelihood of converging to a satisficing solution). Specify fitnessfunc like this: fitnessfunc = @(x) jointhistogramfitness(spmaux, template_vol, input_vol, x);
    % * recombinationrate: crossover rate, the probability that one loser's gene will be overwritten by the winner's.
    % * mutationrate: rate of randomly assigning a random value to one loser's gene.
    % * tournamentcount: number of rounds to find a winner and a loser candidates and update the loser's genes. A high value increases the likelihood of converging to a satisficing solution.
    % * multistart: number of independent rounds to relaunch the whole genetic algorithm, with a brand new population. A high value reduces the variance of the final result. This can be used to test the influence of different parameters values (ie, to hyper-optimize). Multistart allows to reduce variance since we restart from multiple starting points and travel different paths, but this does not help with converging to a better result, it's better to increase tournamentcount than multistartcount to improve performances. Increasing multistart is great to test different parameters (eg, recombinationrate or mutationrate) and see in one result whether the new parameter really improves the performance.
    % * randomstart: defines how the initial population is generated: false will generate a whole population as a copy of the input candidate genestart, true will keep genestart as the reference candidate but all other individuals will be randomly generated according to the randomgen function. If randomstart = true, randomgen needs to be a function accepting popsize as 1st argument and length of genestart as 2nd argument, and which generates popsize individuals to add in the pool. For example, if fitnessfunc = @(x)sum(x), randomgen can be @(popsize,numel_genestart)randi(10, popsize, numel_genestart)
    % * randomgen: function that defines how the initial population is randomly generated when randomstart == true. randomgen should expect 2 parameters: popsize and numel(genestart).
    % * mutatefunc: defines how the mutation is done (ie, how random values are calculated to be assigned to loser's genes). It should expect one parameter which is the locus location i, so that a different mutation scheme can be applied depending on where in the vector we are mutating. If false, will apply a different function (and hence different ranges of values) for translation, rotation and shearing. Example: mutatefunc = @(i)randi(10)
    % * truerandom: defines whether the function is deterministic (false) or non-deterministic (true or a rand('seed') integer)
    % * onlyrigid: if true (default), only rigid-body transforms will be generated (ie, translation and rotation). If false, scale and shearing will also be explored.
    %
    % Outputs:
    % * best: the locally optimal candidate imatrix with all the transform parameters as a 12 items vector (use spm_matrix(best) to transform to an orientation matrix).
    % * bestfit: fitness score of the best candidate.
    % * bestfitstd: standard deviation over all multistart rounds (allows to assess performance variance by evaluating the variability of the current set of parameters).
    % * bestpool and bestpoolfit: all the best candidates across the multistart rounds. This can be used externally to manually select the best candidate using another fitnessfunction (possibly more precise but more time consuming, but since the set is reduced this is faster).
    %
    % Example usage: [best, bestfit, bestfitstd] = mga_imatrix([1 2 3 4 5 6 7 8], 100, 10, @(x)sum(x), 0.7, 0.25, 1000, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
    % another example with a smaller number of tournament rounds and hence smaller population, with similar performances: [best, bestfit, bestfitstd] = mga_imatrix([1 2 3 4 5 6 7 8], 10, 3, @(x)sum(x), 0.7, 0.25, 100, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
    %
    % Tips: The 4 most important parameters to hyper-optimize for accuracy performance are: popsize, demesize, recombinationrate and mutationrate. Assess using (maximizing) bestfit and (minimize) bestfitstd with a high number of (>1000) multistart rounds, which would mean that the genetic algorithm with your parameters and fitnessfunc can reach a high score while having low variance. Also, increasing tournamentcount of course allows the algorithm to converge to a better solution, so optimize your fitnessfunc to be super fast (can use caching for example).
    %
    % Reference: algorithm from: The Microbial Genetic Algorithm, Inman Harvey, 1996
    % There exists a single-line version, see "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534
    %
    % For a simpler implementation, take a look at aux/autogen_mini.m
    %
    % License: MIT License.
    % Copyright (C) 2020-2024 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
    %
    % Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    % The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
        % Default values
        if ~exist('multistart', 'var') || isempty(multistart)
            multistart = 1;
        end
        if ~exist('debugmode', 'var')
            debugmode = false;
        end
        if ~exist('randomstart', 'var') || isempty(randomstart)
            randomstart = false;
        end
        if ~exist('mutatefunc', 'var') || isempty(mutatefunc)
            mutatefunc = false;
        end
        if ~exist('randomgen', 'var') || isempty(randomgen)
            randomgen = false;
        end
        if ~exist('truerandom', 'var') || isempty(randomgen)
            truerandom = false;
        end
        if ~exist('onlyrigid', 'var') || isempty(onlyrigid)
            onlyrigid = true;
        end
    
        % Sanity checks
        if ~isvector(genestart)
            % ensure genestart is a vector
            genestart = spm_imatrix(genestart);
        else
            % else if it's already a vector, convert to a matrix and back to a vector to ensure that all parameters are provided (so that we can simply provide [0] and the rest will be filled)
            genestart = spm_imatrix(spm_matrix(genestart));
        end
    
        % Init multistart gene pool
        bestpool = repmat(genestart, [multistart, 1]);
        bestpoolfit = zeros(multistart, 1);
        if onlyrigid
            genescount = 6;  % restrict to rotation and translation (first 6 items in spm_matrix())
        else
            genescount = numel(genestart);  % explore all genes (rotation, translation, scale and shearing)
        end
        % Memoize fitness function, this will also work between multistart rounds. This was expected to save a lot of time, but in practice it slows down calculations in most cases except in specific circumstances.
        %fitnessfunc_memo = memoize(fitnessfunc);  % deprecated, this slows down a lot the calculations. Use only if your fitness function is computation intensive AND the genes are discrete, else with a continuous valued vector it's useless to memoize
    
        % Save current seed for the random number generator, because SPM functions can reset it (although it is bad practice, their goal is to transform non-deterministic functions to deterministic like that... But this may explain why the same analysis can give different results when ran on different computers/OSes, as this changes the underlying random number generators)
        % We will do 3 things at once here:
        % 1- save the original rng so we can restore it at the end
        orig_rng = rng();
        % 2- set the rng so this function is deterministic
        % Note: setting the generator type is discouraged but we want that so that results are reproducible in the future, we don't need secure pseudorandomness. For more infos, see: https://www.mathworks.com/help/matlab/math/updating-your-random-number-generator-syntax.html
        if isstruct(truerandom)
            %rng(truerandom);
            rand('seed', truerandom);
        elseif ~truerandom
            %rng(0,'philox');  % philox works better than twister. It's important to have a good rng, else the exploration of new solutions by random mutations and random start population will be less effective
            % DEPRECATED: rng() is advised by mathworks but it's magnitude of times slower. It's not made to be in a loop, but we need to reset the rng inside the loop since any call to the fitnessfunc may reset the rng (eg, if fitnessfunc calls spm functions). So we are back to using rand('seed'), which is the fastest option: https://stackoverflow.com/questions/39251206/speed-up-random-number-generation-in-matlab
            % See also these threads about this performance issue of rng():
            % https://www.mathworks.com/matlabcentral/answers/128411-rng-is-slow-to-control-random-number-generation-can-i-go-back-to-seed-or-state
            % https://www.mathworks.com/matlabcentral/answers/5320-speed-improvement-of-the-random-generator
            % https://www.mathworks.com/matlabcentral/answers/67667-performance-degradation-of-random-number-generation-in-matlab-r2013a
            rand('seed', 0);
        % else: we don't reset the rng, we leave as default and hence we get true pseudorandomness (non-deterministic output since we continue with the previous rng state)
        end
        % 3- save the current rng state, so that we can restore it after each spm call
        cur_rng = rand('seed');  %cur_rng = rng();
    
        % For each multistart round (where we restart anew from scratch - this is different from tournaments where each tournament round reuses the same population, here we change the whole population)
        for m=1:multistart
            fprintf('Multistart: %i/%i\n', m, multistart);
            % Init population's genes pool vars
            if ~randomstart
                % Initialize the gene pool by simply replicating the provided genestart vector
                genepool = repmat(genestart, [popsize 1]);
                genepoolfit = ones(popsize,1) .* fitnessfunc(genestart);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
                bestfit = NaN;  % it's useless to compute a bestfit from input because anyway we replicated the input so we are guaranteed to only find a better solution
                best = NaN;
            else
                % Else initialize the gene pool by generating a random translation and orientation for each individual, but retain the original scale and shear from the genestart vector
                if randomgen ~= false
                    % User-defined function
                    genepool = randomgen(popsize, numel(genestart));
                else
                    % Default function
                    if onlyrigid
                        genepool = [round(rand(popsize, 3).*maxdim) - maxdim, rand(popsize,3).*(2*pi) - pi, repmat(genestart(7:end), [popsize 1])];  % random rotation and translation
                        %genepool = [repmat(genestart(1:3), [popsize 1]), rand(popsize,3).*(2*pi) - pi, repmat(genestart(7:end), [popsize 1])];  % random rotation only
                    else
                        genepool = [round(rand(popsize, 3).*maxdim) - maxdim, rand(popsize,3).*(2*pi) - pi, 0.5 + rand(popsize, 6)*1.5];  % random rotation, translation, scale and shearing
                    end
                end  % endif
                genepoolfit = NaN(popsize,1);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
                % but keep the first candidate as the initial one, so we ensure that any candidate we choose is not worse than the input
                genepool(1, :) = genestart;
                best = 1;
                bestfit = fitnessfunc(genestart);
                %bestfit = NaN;
                %best = NaN;
            end  % endif
            % Launch the tournament
            for t=1:tournamentcount
                if mod(t, 10) == 0, fprintf('Tournament: %i/%i\n', t, tournamentcount); end
                % Restore previous state of rng (to avoid SPM meddling with rng - each call to fitnessfunc calls spm and hence meddles with the rng)
                rand('seed', cur_rng);  % rng(cur_rng); -- much slower...
                % Randomly select one individual A
                A = randi(popsize);  % alternative: A = ceil(rand()*popsize);
                % Randomly select another individual B in the deme just after A
                B = mod((A+1+randi(demesize)), popsize)+1;  % alternative: B = mod((A+1+ceil(rand()*demesize)), popsize)+1;
                % Save current rng state (to restore at next tournament)
                cur_rng = rand('seed');  %cur_rng = rng();
                if debugmode, disp([A B]); end;
                % Compute fitness cost for each candidate
                if ~isnan(genepoolfit(A))
                    % If there is a cache for this candidate, use it
                    Afit = genepoolfit(A);
                else
                    % Else compute the fitness cost
                    Afit = fitnessfunc(genepool(A,:));  % memoize fitness to optimize, so that we can reuse directly at the end
                end
                if ~isnan(genepoolfit(B))
                    Bfit = genepoolfit(B);
                else
                    Bfit = fitnessfunc(genepool(B,:));
                end
                % Find the winner
                if debugmode, disp(genepool(A,:)); disp(genepool(B,:)); disp(Afit); disp(Bfit); end
                if (Afit > Bfit)
                    winner = A;
                    loser = B;
                    winnerfit = Afit;
                else
                    winner = B;
                    loser = A;
                    winnerfit = Bfit;
                end  % endif
                % Update fitness cost cache ...
                genepoolfit(loser) = NaN;  % ... by deleting (NaN) the loser's fitness cost, so next time it will be recomputed...
                genepoolfit(winner) = winnerfit;  % ... and by caching the winner's fitness cost.
                % Compare winner with the best fit (memoization)
                if winnerfit >= bestfit || isnan(bestfit)  % note: it's crucial to use >= (and not >) because it can happen that both candidates are equal (they maxxed out), in this case there will still be a loser who will be mutated, which can hence become suboptimal. In that case, if the loser was the best candidate, we need to pass the best label to the winner (despite them being equal), because the winner will not change.
                    bestfit = winnerfit;
                    best = winner;
                    %disp([best, bestfit, genepool(best, :)]);  % debugline
                end  % endif
                % Recombine and mutate for each gene
                rand('seed', cur_rng);  %rng(cur_rng);  % Restore previous state of rng (to avoid SPM meddling with rng)
                for i=1:genescount
                    r = rand();
                    if r < (recombinationrate+mutationrate)  % optimization, see slide 20 of: https://fr.slideshare.net/lrq3000/pathway-evolution-algorithm-in-netlogo
                        if r < recombinationrate
                            % Recombine/crossover (ie, take the allele from the winner)
                            genepool(loser, i) = genepool(winner, i);
                        else
                            % Mutate
                            if mutatefunc ~= false
                                % User-defined mutation function
                                genepool(loser, i) = mutatefunc(i);
                            else
                                % Default mutation function
                                if 1 <= i && i <= 3
                                    % Translation, constrained to image dimension
                                    genepool(loser, i) = round(rand()*(maxdim(i).*2) - maxdim(i));
                                elseif 4 <= i && i <= 6
                                    % Rotation, constrained to radians
                                    genepool(loser, i) = rand()*2*pi - pi;
                                elseif i >= 7
                                    % Scaling and shearing, bounded inside 0.5 to 2.0 here (but in reality it's not bounded, but then the exploration space is way to big)
                                    genepool(loser, i) = 0.5 + rand()*1.5;
                                end  % endif
                            end
                        end  % endif
                    end  % endif
                end  % endfor genes walking
                cur_rng = rand('seed');  %cur_rng = rng(); % Save current rng state (to restore at next tournament)
            end  % endfor tournament rounds
    
            % Select the best candidate
            % DEPRECATED: manual comparison by iterating and comparing all candidates. This does not use memoization. If fitnessfunc is expensive, this will take a long time to compute
            %best_bef = best;
            %bestfit_bef = bestfit;
            %best = 1;
            %bestfit = fitnessfunc(genepool(best, :));
            %for i=2:popsize
            %    newfit = fitnessfunc(genepool(i, :));
            %    if newfit > bestfit
            %        best = i;
            %        bestfit = newfit;
            %    end  % endif
            %end  % endif
            % disp([best_bef best bestfit_bef bestfit]);  % debugline
    
            %disp([best bestfit]);  % debugline
            %disp(genepool(best, :));  % debugline
    
            % Save best candidate of this multistart run
            bestpool(m,:) = genepool(best, :);
            bestpoolfit(m) = bestfit;
        end  %endfor multistart rounds
    
        % End of multistart: select the best candidate over all runs
        [~, bestidx] = max(bestpoolfit);
        best = bestpool(bestidx, :);
        bestfit = bestpoolfit(bestidx);
        bestfitstd = std(bestpoolfit);
    
        % Restore original rng
        rng(orig_rng);
    end  % endfunction
