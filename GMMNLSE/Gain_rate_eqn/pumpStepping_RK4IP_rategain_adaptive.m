function [P1, p5,...
          opt_deltaZ,success] = pumpStepping_RK4IP_rategain_adaptive(sim, gain_rate_eqn,...
                                                                     P0, p5_1,...
                                                                     omegas, dt,...
                                                                     dummy_var)
%PUMPSTEPPING_RK4IP_RATEGAIN_ADAPTIVE Take one step with RK4IP with a gain 
%model solved from rate equations. The gain term is treated as a dispersion
%term, instead of a nonlinear term.
%
% Input:
%    sim.deltaZ - step size; m
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.SK_factor - SK = SK_factor * fiber.SR
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.Raman_model - which Raman model is used
%    sim.Raman_sponRS - consider spontaneous Raman or not
%
%    gain_rate_eqn - container of rate-eqn-gain parameters
%
%    P0 - the counter-propagating pump power
%    p5_1 - the previous RK4 term that can be reused here
%
%    omegas - angular frequencies in 1/ps, in the fft ordering
%    dt - time grid point spacing; ps
%
%    cross_sections_pump
%    cross_sections
%    overlap_factor - no unit for single-mode and 1/um^2 for multimode
%    N_total - (Nx,Nx); the doped ion density; in "1/um^3"
%    FmFnN - the integral2(overlap_factor*N_total) for the signal and ASE
%    GammaN - the integral2(overlap_factor*N_total) for the pump
%
%    dummy_var - unused variable
%
% Output:
%    P1 - the counter-propagating pump power
%    p5 - the RK4 term that can be reused in the next step
%    opt_deltaZ - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

% Propagate through the nonlinearity
if isempty(p5_1)
    p5_1 = N_op(sim, gain_rate_eqn,...
                P0,...
                omegas, dt,...
                dummy_var);
end
p1 = p5_1;
p2 = N_op(sim, gain_rate_eqn,...
          P0+p1*(sim.deltaZ/2),...
          omegas, dt,...
          dummy_var);
p3 = N_op(sim, gain_rate_eqn,...
          P0+p2*(sim.deltaZ/2),...
          omegas, dt,...
          dummy_var);
p4 = N_op(sim, gain_rate_eqn,...
          P0+p3*(sim.deltaZ),...
          omegas, dt,...
          dummy_var);

P1 = P0 + (p1+2*p2+2*p3)*(sim.deltaZ/6) + p4*(sim.deltaZ/6);

% Local error estimate
p5 = N_op(sim, gain_rate_eqn,...
          P1,...
          omegas, dt,...
          dummy_var);
err = sum(abs((p4-p5)*(sim.deltaZ/10)).^2,1);

% Stepsize control
if isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_deltaZ = 0.5*sim.deltaZ;
    success = false;
else
    opt_deltaZ = max(0.5,min(2,0.8*(sim.adaptive_deltaZ.threshold/10/err)^(1/4)))*sim.deltaZ;

    success = err < sim.adaptive_deltaZ.threshold;
end

end

function dPdz = N_op(sim, gain_rate_eqn,...
                          Power_pump_backward,...
                          omegas, dt,...
                          dummy_var)
%N_op Calculate dPdz for the backward pump power

new_Power_pump_backward = solve_gain_rate_eqn('backward',...
                                              sim,gain_rate_eqn,...
                                              dummy_var,dummy_var,0,Power_pump_backward,dummy_var,dummy_var,...
                                              omegas,dt,...
                                              true);
dPdz = (new_Power_pump_backward - Power_pump_backward)/sim.deltaZ;

end