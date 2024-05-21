function gain_rate_eqn = gain_medium(gain_rate_eqn)
%GAIN_MEDIUM Summary of this function goes here
%   Detailed explanation goes here

switch gain_rate_eqn.gain_medium
    case 'Yb'
        gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
        gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
    case 'Er'
        gain_rate_eqn.cross_section_filename = 'Er (optiwave and Smirnov).txt';
        gain_rate_eqn.tau = 9e-3; % lifetime of Er in I_(13/2) state
    case 'Nd'
        gain_rate_eqn.cross_section_filename = 'Nd (Martin thesis 2006).txt';
        gain_rate_eqn.tau = 420e-6; % lifetime of Nd in F_(3/2) state (Bartolacci et al., "Effects of ions clustering in Nd3+/Al3+-codoped double-clad fiber laser operating near 930 nm")
        gain_rate_eqn.cluster_fraction = 0.2;
    case 'Tm'
        
    case 'Ho'
        
end

end

