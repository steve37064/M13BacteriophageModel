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
%   species_init    : row vector of 72 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 64 model parameters.
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
   parameters = [ 0.1, 0.117, 1.0, 1.0E2, 1, 1E-15, 1E-3, 1100, 1.0, 6E-3, 1E-2, 6.5E-2, 2E-3, 2E-4, 1.34, 3.6E-2, 7.2E-2, 2.5E-1, 5.5E-2, 4.3E-2, 0.60, 2.0E-2, 5.8E-2, 4.6E-3, 1.9E-3, 1.4E-3, 1.2E-3, 1.9E-3, 1.9E-3, 1.9E-3, 1.60, 10.8, 8.00, 4.60, 20.0, 0.80, 4.40, 1.60, 1.08, 4.00, 0.25, 4.8E-2, 4.1E-2, 3.7E-2, 3.8E-2, 4.3E-1, 2.5E-1, 8.9E-1, 2.4E-1, 3.1E-1, 80, 1.5E-4, 4, 1.0E-2, 2.6E-5, 1.0, 0.1, 0.1, 0.1, 1.2E-3, 0.1, 3000, 40, 1E-6 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 64  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 64].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 72  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 72].\n' );
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
param_labels = { 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'Km1', 'n1', 'C8_A', 'C8_B', 'C8_H', 'C8_Z', 'C8_W', 'C9', 'C10_A', 'C10_B', 'C10_H', 'C10_Z', 'C10_W', 'C11', 'C12_A', 'C12_D', 'C12_E', 'C12_F', 'C12_G', 'C12_H', 'C12_W', 'C12_Y', 'C12_Z', 'C13_1', 'C13_2', 'C13_3', 'C13_4', 'C13_5', 'C13_6', 'C13_8', 'C13_9', 'C13_10', 'C13_11', 'C14', 'C15_1', 'C15_2', 'C15_3', 'C15_4', 'C15_5', 'C15_6', 'C15_8', 'C15_10', 'C15_11', 'EF1', 'EF2', 'EF3', 'EF5', 'EF10', 'n2', 'C16', 'C17', 'C18', 'C19', 'C20', 'Km2', 'n3', 'EPS' };



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
observables_out = zeros( length(timepoints), 72 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'DP3', 'RNAP', 'R', 'ssDNA', 'ssPDNA', 'RF1', 'RF2', 'RF2DP3', 'DA', 'DB', 'DH', 'DZ', 'DW', 'EA', 'EB', 'EH', 'EZ', 'EW', 'ELA', 'ELB', 'ELH', 'ELZ', 'ELW', 'A', 'D', 'E', 'F', 'G', 'H', 'W', 'Y', 'Z', 'RBS2', 'RBS5', 'RBS9', 'RBS8', 'RBS1', 'RBS3', 'RBS4', 'RBS6', 'RBS11', 'RBS1R', 'RBS2R', 'RBS3R', 'RBS4R', 'RBS5R', 'RBS6R', 'RBS8R', 'RBS9R', 'RBS10R', 'RBS11R', 'PD1', 'PD2', 'PD3', 'PD4', 'PD5', 'PD6', 'PD8', 'PD9', 'PD10', 'PD11', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11' };

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

    species_init = zeros(1,72);
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
    species_init(29) = 0;
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;
    species_init(54) = 0;
    species_init(55) = 0;
    species_init(56) = 0;
    species_init(57) = 0;
    species_init(58) = 0;
    species_init(59) = 0;
    species_init(60) = 0;
    species_init(61) = 0;
    species_init(62) = 0;
    species_init(63) = 0;
    species_init(64) = 0;
    species_init(65) = 0;
    species_init(66) = 0;
    species_init(67) = 0;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 0;
    species_init(71) = 0;
    species_init(72) = 0;

end


% user-defined functions
% function P5InhibitionP2
function [val] = P5InhibitionP2(expressions, observables)
    val = (1/(1+((observables(66)/(expressions(52)*(observables(33)+1)))^expressions(56))));
end

% function RBS2Removal
function [val] = RBS2Removal(expressions, observables)
    val = (1/(observables(33)+expressions(64)));
end

% function RBS5Removal
function [val] = RBS5Removal(expressions, observables)
    val = (1/(observables(34)+expressions(64)));
end

% function RBS9Removal
function [val] = RBS9Removal(expressions, observables)
    val = (1/(observables(35)+expressions(64)));
end

% function RBS8Removal
function [val] = RBS8Removal(expressions, observables)
    val = (1/(observables(36)+expressions(64)));
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

    expressions = zeros(1,67);
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
    expressions(65) = (0.7*expressions(23));
    expressions(66) = (0.7*expressions(24));
    expressions(67) = (0.7*expressions(26));
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,72);
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
    observables(29) = species(29);
    observables(30) = species(30);
    observables(31) = species(31);
    observables(32) = species(32);
    observables(33) = species(33);
    observables(34) = species(34);
    observables(35) = species(35);
    observables(36) = species(36);
    observables(37) = species(37);
    observables(38) = species(38);
    observables(39) = species(39);
    observables(40) = species(40);
    observables(41) = species(41);
    observables(42) = species(42);
    observables(43) = species(43);
    observables(44) = species(44);
    observables(45) = species(45);
    observables(46) = species(46);
    observables(47) = species(47);
    observables(48) = species(48);
    observables(49) = species(49);
    observables(50) = species(50);
    observables(51) = species(51);
    observables(52) = species(52);
    observables(53) = species(53);
    observables(54) = species(54);
    observables(55) = species(55);
    observables(56) = species(56);
    observables(57) = species(57);
    observables(58) = species(58);
    observables(59) = species(59);
    observables(60) = species(60);
    observables(61) = species(61);
    observables(62) = species(62);
    observables(63) = species(63);
    observables(64) = species(64);
    observables(65) = species(65);
    observables(66) = species(66);
    observables(67) = species(67);
    observables(68) = species(68);
    observables(69) = species(69);
    observables(70) = species(70);
    observables(71) = species(71);
    observables(72) = species(72);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,72);
    ratelaws(1) = expressions(1)*species(4)*species(1);
    ratelaws(2) = expressions(2)*species(5);
    ratelaws(3) = expressions(3)*species(6)*species(63);
    ratelaws(4) = expressions(4)*species(7);
    ratelaws(5) = rateLaw__1(expressions,observables)*species(7)*species(1)*species(9)*species(10)*species(11)*species(12)*species(13);
    ratelaws(6) = expressions(2)*species(8);
    ratelaws(7) = expressions(10)*species(9)*species(2);
    ratelaws(8) = expressions(15)*species(14);
    ratelaws(9) = expressions(16)*species(19);
    ratelaws(10) = rateLaw__2(expressions,observables)*species(24)*species(33);
    ratelaws(11) = (0.7*expressions(23))*species(25);
    ratelaws(12) = rateLaw__4(expressions,observables)*species(25)*species(34)*species(35)*species(36);
    ratelaws(13) = (0.7*expressions(24))*species(26);
    ratelaws(14) = rateLaw__6(expressions,observables)*species(26)*species(34)*species(35)*species(36);
    ratelaws(15) = rateLaw__7(expressions,observables)*species(27)*species(34);
    ratelaws(16) = rateLaw__8(expressions,observables)*species(27)*species(34)*species(35)*species(36);
    ratelaws(17) = (0.7*expressions(26))*species(28);
    ratelaws(18) = rateLaw__10(expressions,observables)*species(28)*species(35)*species(36);
    ratelaws(19) = rateLaw__11(expressions,observables)*species(29)*species(35)*species(36);
    ratelaws(20) = rateLaw__12(expressions,observables)*species(33)*species(3);
    ratelaws(21) = expressions(41)*species(43);
    ratelaws(22) = expressions(43)*species(53);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(72,1);
    
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
    Dspecies(15) = 0.0;
    Dspecies(16) = 0.0;
    Dspecies(17) = 0.0;
    Dspecies(18) = 0.0;
    Dspecies(19) = ratelaws(8) -ratelaws(9);
    Dspecies(20) = 0.0;
    Dspecies(21) = 0.0;
    Dspecies(22) = 0.0;
    Dspecies(23) = 0.0;
    Dspecies(24) = ratelaws(9) -ratelaws(10);
    Dspecies(25) = ratelaws(10) -ratelaws(11) -ratelaws(12);
    Dspecies(26) = ratelaws(11) -ratelaws(13) -ratelaws(14);
    Dspecies(27) = ratelaws(13) -ratelaws(15) -ratelaws(16);
    Dspecies(28) = ratelaws(15) -ratelaws(17) -ratelaws(18);
    Dspecies(29) = ratelaws(17) -ratelaws(19);
    Dspecies(30) = 0.0;
    Dspecies(31) = 0.0;
    Dspecies(32) = 0.0;
    Dspecies(33) = ratelaws(9) -ratelaws(10) -ratelaws(20) +ratelaws(21);
    Dspecies(34) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(15) -ratelaws(16);
    Dspecies(35) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(16) -ratelaws(18) -ratelaws(19);
    Dspecies(36) = ratelaws(9) -ratelaws(12) -ratelaws(14) -ratelaws(16) -ratelaws(18) -ratelaws(19);
    Dspecies(37) = 0.0;
    Dspecies(38) = 0.0;
    Dspecies(39) = 0.0;
    Dspecies(40) = 0.0;
    Dspecies(41) = 0.0;
    Dspecies(42) = 0.0;
    Dspecies(43) = ratelaws(20) -ratelaws(21);
    Dspecies(44) = 0.0;
    Dspecies(45) = 0.0;
    Dspecies(46) = 0.0;
    Dspecies(47) = 0.0;
    Dspecies(48) = 0.0;
    Dspecies(49) = 0.0;
    Dspecies(50) = 0.0;
    Dspecies(51) = 0.0;
    Dspecies(52) = 0.0;
    Dspecies(53) = ratelaws(21) -ratelaws(22);
    Dspecies(54) = 0.0;
    Dspecies(55) = 0.0;
    Dspecies(56) = 0.0;
    Dspecies(57) = 0.0;
    Dspecies(58) = 0.0;
    Dspecies(59) = 0.0;
    Dspecies(60) = 0.0;
    Dspecies(61) = 0.0;
    Dspecies(62) = 0.0;
    Dspecies(63) = -ratelaws(3) +ratelaws(4) +ratelaws(6) +ratelaws(22);
    Dspecies(64) = 0.0;
    Dspecies(65) = 0.0;
    Dspecies(66) = 0.0;
    Dspecies(67) = 0.0;
    Dspecies(68) = 0.0;
    Dspecies(69) = 0.0;
    Dspecies(70) = 0.0;
    Dspecies(71) = 0.0;
    Dspecies(72) = 0.0;

end


end
