function [err, timepoints, species_out, observables_out ] = Model( timepoints, species_init, parameters, suppress_plot )
%MODEL Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Model' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file Model_cvode.c must be compiled (see that file for details).
%   MODEL returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Model( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of 81 initial species populations.
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
   parameters = [ 0.1, 0.125, 1.0, 1.0E3, 1, 9.5454e-16, 1E-3, 1100, 1.0, 6E-3, 1E-2, 6.5E-3, 2E-4, 2.25E-3, 1.34, 0.036216, 0.070526, 0.24815, 0.054918, 0.043226, 0.60, 0.020387, 0.0057762, 0.0057763, 0.004621, 0.00192542, 0.00144402, 0.0011552, 0.0019254, 0.0019254, 0.0019254, 1.60, 10.8, 8.00, 4.60, 20.0, 0.80, 4.40, 1.60, 1.08, 4.00, 0.25, 0.048443, 0.040896, 0.037267, 0.038286, 0.43299, 0.2515, 0.89362, 0.23729, 0.30657, 80, 1.5E-4, 4, 1.0E-2, 2.6E-5, 1.0, 0.1, 0.1, 0.1, 1.23E-3, 0.1, 3000, 40, 1E-6 ];
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
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 81  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 81].\n' );
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
param_labels = { 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'Km1', 'n1', 'C8_A', 'C8_B', 'C8_H', 'C8_Z', 'C8_W', 'C9', 'C10_A', 'C10_B', 'C10_H', 'C10_W', 'C10_Z', 'C11', 'C12_A', 'C12_B', 'C12_D', 'C12_E', 'C12_F', 'C12_G', 'C12_H', 'C12_W', 'C12_Y', 'C12_Z', 'C13_1', 'C13_2', 'C13_3', 'C13_4', 'C13_5', 'C13_6', 'C13_8', 'C13_9', 'C13_10', 'C13_11', 'C14', 'C15_1', 'C15_2', 'C15_3', 'C15_4', 'C15_5', 'C15_6', 'C15_8', 'C15_10', 'C15_11', 'EF1', 'EF2', 'EF3', 'EF5', 'EF10', 'n2', 'C16', 'C17', 'C18', 'C19', 'C20', 'Km2', 'n3', 'EPS' };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = Model_cvode( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'DP3', 'RNAP', 'R', 'ssDNA', 'ssPDNA', 'RF1', 'RF2', 'RF2DP3', 'P5DNA', 'DA', 'DB', 'DH', 'DZ', 'DW', 'EA', 'EB', 'EH', 'EZ', 'EW', 'ELA', 'ELB', 'ELH', 'ELZ', 'ELW', 'A', 'B', 'D', 'E', 'F', 'G', 'H', 'W', 'Y', 'Z', 'RBS2', 'RBS10', 'RBS5', 'RBS9', 'RBS8', 'RBS1', 'RBS3', 'RBS4', 'RBS6', 'RBS11', 'RBS1R', 'RBS2R', 'RBS3R', 'RBS4R', 'RBS5R', 'RBS6R', 'RBS8R', 'RBS9R', 'RBS10R', 'RBS11R', 'PD1', 'PD2', 'PD3', 'PD4', 'PD5', 'PD6', 'PD8', 'PD9', 'PD10', 'PD11', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P2P10', 'As', 'PI', 'PE', 'PF', 'Phage' };

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



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,81);
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
    species_init(67) = 5;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 5;
    species_init(71) = 5;
    species_init(72) = 2700;
    species_init(73) = 5;
    species_init(74) = 0;
    species_init(75) = 0;
    species_init(76) = 0;
    species_init(77) = 0;
    species_init(78) = 0;
    species_init(79) = 0;
    species_init(80) = 0;
    species_init(81) = 0;

end


end
