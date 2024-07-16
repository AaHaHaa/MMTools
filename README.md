# MMTools
This is the shared package to simulate pulse propagation in a solid-core fiber with GMMNLSE/MM-UPPE with MATLAB.

It is useful for simulating single-mode/multimode mode-locking/oscillators, fiber amplifiers, single-mode/vector/multimode solitons, spatial beam cleaning in multimode fibers, fiber optical parametric amplifier (FOPA), and so. Some typical examples of oscillators include all-normal-dispersion (ANDi) oscillators and Mamyshev oscillators. Amplifiers include linear chirped-pulse amplification (CPA) and gain-managed nonlinear amplification (GMNA).

## Capabilities:<br>
1. It solves the pulse propagation with
   - RK4IP if single mode
   - MPA if multimode
2. Support both scalar and polarized scenarios. Controlled with `sim.scalar=true/false`.
3. Adaptive step-size control are implemented for both RK4IP and MPA.

> Note:<br>
Although adaptive-step-size control for RK4IP isn't new with published papers, adaptive-step-size control for MPA is new. I didn't publish a separate paper discussing this numerical scheme, which is perhaps the fastest and the most convenient numerical scheme for general multimode situations (as in step-index, graded-index, or hollow-core fibers, etc.) by far (written on 2/14/2024). The implementation detail is described in the supplement of https://doi.org/10.1364/JOSAB.500586.

4. Support random mode coupling.
5. Support both passive and gain fibers
   - Gain model includes Gaussian gain and rate-equation gain, for both single-mode and multimode scenarios.
   - For rate-equation-gain modeling, all pumping schemes are implemented: copumping, counterpumping, co+counter-pumping, as well as with and without ASE.
   - If ASE is included, the effect of ASE to the coherent signal field is simulated, rather than only a separate power variable $P_{\text{ASE}}(\omega)$ from the coherent signal field ${F[A]}(\omega)$.
   - Rate-equation model supports Nd, Yb, Er, Tm, Ho. For more details, see `readme.pdf`.
   - Support ring- and linear-oscillator configurations with fast convergence (with the use of `saved_data`). For linear oscillators, inclusion of influence from pulses of both directions to the gain medium is considered. As an example, please see the numerical section of http://josab.osa.org/abstract.cfm?URI=josab-38-3-743 to understand the necessity of this two-pulse saturation effect in a linear oscillator.
6. Support the addition of spontaneous Raman scattering and input-pulse shot noise. For details, see the supplement of https://doi.org/10.1063/5.0189749. Although this paper is for gas-filled HCF, its Raman formulation actually works generally anywhere, such as solid-core fibers. It's just a mathematical description of how to treat the sinusoidal Raman response under different temporal scales.
7. For multimode, GPU computations (with Nvidia CUDA) is highly recommended. I have written a lot of CUDA files to speed up simulations. It is controlled by `sim.gpu_yes=true/false`.

## Notes:<br>
For details, please read the supplement of our paper: https://doi.org/10.1364/JOSAB.500586.  
Please don't forget to cite our paper if you find this code useful in your work. I, the young and early-career researcher, need your support. Similarly, if you need help or have questions about the code, please feel free to send me an email.

There is a `readme.pdf` in the `Documentations/` folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics. This code is basically an upgraded and highly-optimized version of our "WiseLabAEP/GMMNLSE-Solver-FINAL" from "https://github.com/WiseLabAEP/GMMNLSE-Solver-FINAL," with much more functionalities, which however might overwhelm users and thus require more fiber-optic background. It can run order-of-magnitude faster than our old code due to optimizing with CUDA+shared memory, as well as reducing the usage of for-loops. Although our old one claims to be fast with GPU, its CUDA implementation is not optimized, let alone its CPU implementation with a lot of slow for-loops. Besides, this package includes adaptive step-size control, which improves the performance significantly and allows users to be free from worrying the reliability of a simulation. For optimization details, please see the supplement of our paper mentioned previously. 

If you have questions, feel free to ask them here or send me an email (email address is in my paper).

## History:<br>
* 11/1/2023:<br>
If you downloaded the code earlier, please re-download it. There was a huge bug in polarization modes. I fixed it only recently. Now it works correctly.
* 12/22/2023:<br>
I've been constantly modifying it, especially the ASE part, these days.<br>
If you find a bug, download the latest one since I've been changing it often. If it still exists, let me know.
* 1/3/2024:<br>
I've fixed the bugs regarding spontaneous Raman scattering. My previous implementation was wrong and created almost negligible spontaneous Raman. It didn't affect any of the results that didn't care about spontaneous Raman scattering.
About the implementation details, please find them in https://doi.org/10.1063/5.0189749.
* 1/17/2024:<br>
Since I've received questions about the Fourier Transform, I've added explanation about it in the readme.pdf. Because of the laser-field definition, Fourier Transform should be `ifft` in MATLAB; be careful about this! It's different from the mathematical convention. This affects phase results and even critical (and can make the result wrong) Fourier-Transform-constant issues, such as different constants of convolution theorem for different conventions.
* 2/14/2024:<br>
There is a significant bug in CUDA related to spontaneous Raman scattering that will simply fail MATLAB. I fixed it. Please download the new one if you downloaded the code between 1/3/2024, when I claimed to add and fix the spontaneous Raman scattering, and today.
* 5/20/2024:<br>
Finish the initial implementation of Er and Nd rate-equation gain modeling. More tests will be done for verification.
* 7/17/2024:<br>
I've fixed bugs related to multimode mode-locking. Thanks Yi Zhou, from Univeristy of Hong Kong, for asking me add examples for a few multimode functions. Please check the "MM ANDi" example in "ANDi oscillator/" folder in "Examples/". In addition, I've finished implementing all types of gain media. Please take a look. More tests need to be done.  
Addition of ASE to the coherent signal field is corrected, which was wrong previously. See the comments in the `stepping_RK4IP/MPA_rategain.m` for details.
