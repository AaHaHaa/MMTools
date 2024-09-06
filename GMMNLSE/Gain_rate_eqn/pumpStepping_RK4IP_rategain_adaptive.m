function [P1, p5,...
          N,...
          opt_dz,success] = pumpStepping_RK4IP_rategain_adaptive(sim, gain_rate_eqn,...
                                                                 P0, p5_1,...
                                                                 N0,...
                                                                 omegas, dt,...
                                                                 dummy_var)
%PUMPSTEPPING_RK4IP_RATEGAIN_ADAPTIVE Take one step with RK4IP with a gain 
%model solved from rate equations. The gain term is treated as a dispersion
%term, instead of a nonlinear term.
%
% Input:
%    sim.dz - step size; m
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.SK_factor - SK = SK_factor * fiber.SR
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.include_Raman - whether to consider Raman
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
%    opt_dz - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance
%
% For adaptive-step implementation, check http://www.sciencedirect.com/science/article/pii/S0010465512004262
%
%    Balac and Mahe, Embedded Runge-Kutta scheme for step-size control in 
%    the interaction picture method, "Comput. Phys. Commun. 184(4), 1211-
%    1219 (2013)

% Propagate through the nonlinearity
if isempty(p5_1)
    p5_1 = N_op(sim, gain_rate_eqn,...
                P0,N0,...
                omegas, dt,...
                dummy_var);
end
p1 = p5_1;
p2 = N_op(sim, gain_rate_eqn,...
          P0+p1*(sim.dz/2),N0,...
          omegas, dt,...
          dummy_var);
p3 = N_op(sim, gain_rate_eqn,...
          P0+p2*(sim.dz/2),N0,...
          omegas, dt,...
          dummy_var);
p4 = N_op(sim, gain_rate_eqn,...
          P0+p3*(sim.dz),N0,...
          omegas, dt,...
          dummy_var);

P1 = P0 + (p1+2*p2+2*p3)*(sim.dz/6) + p4*(sim.dz/6);

% Local error estimate
[p5,N] = N_op(sim, gain_rate_eqn,...
              P1,N0,...
              omegas, dt,...
              dummy_var);
err = sum(abs((p4-p5)*(sim.dz/10)).^2,1);

% Stepsize control
if isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/10/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

end

function [dPdz,N] = N_op(sim, gain_rate_eqn,...
                         Power_pump_backward,N0,...
                         omegas, dt,...
                         dummy_var)
%N_op Calculate dPdz for the backward pump power

[new_Power_pump_backward,~,N] = solve_gain_rate_eqn('backward',...
                                                    sim,gain_rate_eqn,...
                                                    N0,...
                                                    dummy_var,dummy_var,0,Power_pump_backward,dummy_var,dummy_var,...
                                                    omegas,dt,...
                                                    true);
dPdz = (new_Power_pump_backward - Power_pump_backward)/sim.dz;

end