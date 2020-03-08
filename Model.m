function [err, timepoints, species_out, observables_out] = Model( timepoints, species_init, parameters, suppress_plot )
%MODEL Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Model' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. MODEL returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Model( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 28 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 65 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 0.1, 0.117, 1.0, 1.0E2, 1, 1E-15, 1E-3, 1100, 1.0, 6E-3, 1E-2, 6.5E-2, 2E-3, 2E-4, 1.34, 3.6E-2, 7.2E-2, 2.5E-1, 5.5E-2, 4.3E-2, 0.60, 2.0E-2, 5.8E-2, 4.6E-3, 1.9E-3, 1.4E-3, 1.2E-3, 1.9E-3, 1.9E-3, 1.9E-3, 1.60, 10.8, 8.00, 4.60, 20.0, 0.80, 4.40, 1.60, 1.08, 4.00, 0.25, 4.8E-2, 4.1E-2, 3.7E-2, 3.8E-2, 4.3E-1, 2.5E-1, 8.9E-1, 2.4E-1, 3.1E-1, 80, 1.5E-4, 4, 1.0E-2, 2.6E-5, 1.0, 0.1, 0.1, 0.1, 1.2E-3, 0.1, 3000, 40, 1E-6, 0 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 65  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 65].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 28  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 28].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'Km1', 'n1', 'C8_A', 'C8_B', 'C8_H', 'C8_Z', 'C8_W', 'C9', 'C10_A', 'C10_B', 'C10_H', 'C10_Z', 'C10_W', 'C11', 'C12_A', 'C12_D', 'C12_E', 'C12_F', 'C12_G', 'C12_H', 'C12_W', 'C12_Y', 'C12_Z', 'C13_1', 'C13_2', 'C13_3', 'C13_4', 'C13_5', 'C13_6', 'C13_8', 'C13_9', 'C13_10', 'C13_11', 'C14', 'C15_1', 'C15_2', 'C15_3', 'C15_4', 'C15_5', 'C15_6', 'C15_8', 'C15_10', 'C15_11', 'EF1', 'EF2', 'EF3', 'EF5', 'EF10', 'n2', 'C16', 'C17', 'C18', 'C19', 'C20', 'Km2', 'n3', 'EPS', 'P5' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   0.0001,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 28 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'DP3', 'RNAP', 'R', 'ssDNA', 'ssPDNA', 'RF1', 'RF2', 'RF2DP3', 'DA', 'DB', 'DH', 'DZ', 'DW', 'EA', 'ELA', 'A', 'D', 'E', 'F', 'G', 'H', 'RBS2', 'RBS5', 'RBS9', 'RBS8', 'RBS2R', 'PD2', 'P2' };

    % construct figure
    plot(timepoints,observables_out);
    title('Model observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,28);
    species_init(1) = 3.0;
    species_init(2) = 1280.0;
    species_init(3) = 7880.0;
    species_init(4) = 1.0;
    species_init(5) = 0;
    species_init(6) = 0;
    species_init(7) = 0;
    species_init(8) = 0;
    species_init(9) = 0;
    species_init(10) = 0;
    species_init(11) = 0;
    species_init(12) = 0;
    species_init(13) = 0;
    species_init(14) = 0;
    species_init(15) = 0;
    species_init(16) = 0;
    species_init(17) = 0;
    species_init(18) = 0;
    species_init(19) = 0;
    species_init(20) = 0;
    species_init(21) = 0;
    species_init(22) = 0;
    species_init(23) = 0;
    species_init(24) = 0;
    species_init(25) = 0;
    species_init(26) = 0;
    species_init(27) = 0;
    species_init(28) = 0;

end


% user-defined functions
% function P5InhibitionP2
function [val] = P5InhibitionP2(expressions, observables)
    val = (1/(1+((expressions(65)/(expressions(52)*(observables(22)+1)))^expressions(56))));
end

% function RBS2Removal
function [val] = RBS2Removal(expressions, observables)
    val = (1/(observables(22)+expressions(64)));
end

% function RBS5Removal
function [val] = RBS5Removal(expressions, observables)
    val = (1/(observables(23)+expressions(64)));
end

% function RBS9Removal
function [val] = RBS9Removal(expressions, observables)
    val = (1/(observables(24)+expressions(64)));
end

% function RBS8Removal
function [val] = RBS8Removal(expressions, observables)
    val = (1/(observables(25)+expressions(64)));
end

% function DARemoval
function [val] = DARemoval(expressions, observables)
    val = (1/(observables(9)+expressions(64)));
end

% function DBRemoval
function [val] = DBRemoval(expressions, observables)
    val = (1/(observables(10)+expressions(64)));
end

% function DHRemoval
function [val] = DHRemoval(expressions, observables)
    val = (1/(observables(11)+expressions(64)));
end

% function DZRemoval
function [val] = DZRemoval(expressions, observables)
    val = (1/(observables(12)+expressions(64)));
end

% function DWRemoval
function [val] = DWRemoval(expressions, observables)
    val = (1/(observables(13)+expressions(64)));
end

% function rateLaw__1
function [val] = rateLaw__1(expressions, observables)
    val = (((((expressions(1)*DARemoval(expressions,observables))*DBRemoval(expressions,observables))*DHRemoval(expressions,observables))*DZRemoval(expressions,observables))*DWRemoval(expressions,observables));
end

% function rateLaw__2
function [val] = rateLaw__2(expressions, observables)
    val = (expressions(22)*RBS2Removal(expressions,observables));
end

% function rateLaw__4
function [val] = rateLaw__4(expressions, observables)
    val = ((((0.3*expressions(23))*RBS5Removal(expressions,observables))*RBS9Removal(expressions,observables))*RBS8Removal(expressions,observables));
end

% function rateLaw__6
function [val] = rateLaw__6(expressions, observables)
    val = ((((0.3*expressions(24))*RBS5Removal(expressions,observables))*RBS9Removal(expressions,observables))*RBS8Removal(expressions,observables));
end

% function rateLaw__7
function [val] = rateLaw__7(expressions, observables)
    val = ((0.7*expressions(25))*RBS5Removal(expressions,observables));
end

% function rateLaw__8
function [val] = rateLaw__8(expressions, observables)
    val = ((((0.3*expressions(25))*RBS5Removal(expressions,observables))*RBS9Removal(expressions,observables))*RBS8Removal(expressions,observables));
end

% function rateLaw__10
function [val] = rateLaw__10(expressions, observables)
    val = (((0.3*expressions(26))*RBS9Removal(expressions,observables))*RBS8Removal(expressions,observables));
end

% function rateLaw__11
function [val] = rateLaw__11(expressions, observables)
    val = ((expressions(27)*RBS9Removal(expressions,observables))*RBS8Removal(expressions,observables));
end

% function rateLaw__12
function [val] = rateLaw__12(expressions, observables)
    val = (expressions(32)*P5InhibitionP2(expressions,observables));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,68);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = parameters(22);
    expressions(23) = parameters(23);
    expressions(24) = parameters(24);
    expressions(25) = parameters(25);
    expressions(26) = parameters(26);
    expressions(27) = parameters(27);
    expressions(28) = parameters(28);
    expressions(29) = parameters(29);
    expressions(30) = parameters(30);
    expressions(31) = parameters(31);
    expressions(32) = parameters(32);
    expressions(33) = parameters(33);
    expressions(34) = parameters(34);
    expressions(35) = parameters(35);
    expressions(36) = parameters(36);
    expressions(37) = parameters(37);
    expressions(38) = parameters(38);
    expressions(39) = parameters(39);
    expressions(40) = parameters(40);
    expressions(41) = parameters(41);
    expressions(42) = parameters(42);
    expressions(43) = parameters(43);
    expressions(44) = parameters(44);
    expressions(45) = parameters(45);
    expressions(46) = parameters(46);
    expressions(47) = parameters(47);
    expressions(48) = parameters(48);
    expressions(49) = parameters(49);
    expressions(50) = parameters(50);
    expressions(51) = parameters(51);
    expressions(52) = parameters(52);
    expressions(53) = parameters(53);
    expressions(54) = parameters(54);
    expressions(55) = parameters(55);
    expressions(56) = parameters(56);
    expressions(57) = parameters(57);
    expressions(58) = parameters(58);
    expressions(59) = parameters(59);
    expressions(60) = parameters(60);
    expressions(61) = parameters(61);
    expressions(62) = parameters(62);
    expressions(63) = parameters(63);
    expressions(64) = parameters(64);
    expressions(65) = parameters(65);
    expressions(66) = (0.7*expressions(23));
    expressions(67) = (0.7*expressions(24));
    expressions(68) = (0.7*expressions(26));
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,28);
    observables(1) = species(1);
    observables(2) = species(2);
    observables(3) = species(3);
    observables(4) = species(4);
    observables(5) = species(5);
    observables(6) = species(6);
    observables(7) = species(7);
    observables(8) = species(8);
    observables(9) = species(9);
    observables(10) = species(10);
    observables(11) = species(11);
    observables(12) = species(12);
    observables(13) = species(13);
    observables(14) = species(14);
    observables(15) = species(15);
    observables(16) = species(16);
    observables(17) = species(17);
    observables(18) = species(18);
    observables(19) = species(19);
    observables(20) = species(20);
    observables(21) = species(21);
    observables(22) = species(22);
    observables(23) = species(23);
    observables(24) = species(24);
    observables(25) = species(25);
    observables(26) = species(26);
    observables(27) = species(27);
    observables(28) = species(28);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,28);
    ratelaws(1) = expressions(1)*species(4)*species(1);
    ratelaws(2) = expressions(2)*species(5);
    ratelaws(3) = expressions(3)*species(6)*species(28);
    ratelaws(4) = expressions(4)*species(7);
    ratelaws(5) = rateLaw__1(expressions,observables)*species(7)*species(1)*species(9)*species(10)*species(11)*species(12)*species(13);
    ratelaws(6) = expressions(2)*species(8);
    ratelaws(7) = expressions(10)*species(9)*species(2);
    ratelaws(8) = expressions(15)*species(14);
    ratelaws(9) = expressions(16)*species(15);
    ratelaws(10) = rateLaw__2(expressions,observables)*species(16)*species(22);
    ratelaws(11) = (0.7*expressions(23))*species(17);
    ratelaws(12) = rateLaw__4(expressions,observables)*species(17)*species(23)*species(24)*species(25);
    ratelaws(13) = (0.7*expressions(24))*species(18);
    ratelaws(14) = rateLaw__6(expressions,observables)*species(18)*species(23)*species(24)*species(25);
    ratelaws(15) = rateLaw__7(expressions,observables)*species(19)*species(23);
    ratelaws(16) = rateLaw__8(expressions,observables)*species(19)*species(23)*species(24)*species(25);
    ratelaws(17) = (0.7*expressions(26))*species(20);
    ratelaws(18) = rateLaw__10(expressions,observables)*species(20)*species(24)*species(25);
    ratelaws(19) = rateLaw__11(expressions,observables)*species(21)*species(24)*species(25);
    ratelaws(20) = rateLaw__12(expressions,observables)*species(22)*species(3);
    ratelaws(21) = expressions(41)*species(26);
    ratelaws(22) = expressions(43)*species(27);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(28,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1) +ratelaws(2) -ratelaws(5) +ratelaws(6);
    Dspecies(2) = -ratelaws(7) +ratelaws(9);
    Dspecies(3) = -ratelaws(20) +ratelaws(22);
    Dspecies(4) = -ratelaws(1) +ratelaws(6);
    Dspecies(5) = ratelaws(1) -ratelaws(2);
    Dspecies(6) = ratelaws(2) -ratelaws(3) +ratelaws(4) +ratelaws(6);
    Dspecies(7) = ratelaws(3) -ratelaws(4) -ratelaws(5);
    Dspecies(8) = ratelaws(5) -ratelaws(6);
    Dspecies(9) = ratelaws(2) -ratelaws(5) +ratelaws(6) -ratelaws(7) +ratelaws(8);
    Dspecies(10) = ratelaws(2) -ratelaws(5) +ratelaws(6);
    Dspecies(11) = ratelaws(2) -ratelaws(5) +ratelaws(6);
    Dspecies(12) = ratelaws(2) -ratelaws(5) +ratelaws(6);
    Dspecies(13) = ratelaws(2) -ratelaws(5) +ratelaws(6);
    Dspecies(14) = ratelaws(7) -ratelaws(8);
    Dspecies(15) = ratelaws(8) -ratelaws(9);
    Dspecies(16) = ratelaws(9) -ratelaws(10);
    Dspecies(17) = ratelaws(10) -ratelaws(11) -ratelaws(12);
    Dspecies(18) = ratelaws(11) -ratelaws(13) -ratelaws(14);
    Dspecies(19) = ratelaws(13) -ratelaws(15) -ratelaws(16);
    Dspecies(20) = ratelaws(15) -ratelaws(17) -ratelaws(18);
    Dspecies(21) = ratelaws(17) -ratelaws(19);
    Dspecies(22) = ratelaws(9) -ratelaws(10) -ratelaws(20) +ratelaws(21);
    Dspecies(23) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(15) -ratelaws(16);
    Dspecies(24) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(16) -ratelaws(18) -ratelaws(19);
    Dspecies(25) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(16) -ratelaws(18) -ratelaws(19);
    Dspecies(26) = ratelaws(20) -ratelaws(21);
    Dspecies(27) = ratelaws(21) -ratelaws(22);
    Dspecies(28) = -ratelaws(3) +ratelaws(4) +ratelaws(6) +ratelaws(22);

end


end
