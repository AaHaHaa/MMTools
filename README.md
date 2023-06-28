# MMTools
This is the shared package to simulate pulse propagation in a fiber with GMMNLSE/UPPE.

It solves the pulse propagation with RK4IP if single mode and MPA if multimode. Besides, it is implemented with an adaptive step-size control for both methods. Both passive and gain fibers can be simulated, where gain model includes Gaussian gain and rate-equation gain, for both single mode and multimode. Random mode coupling is included as well.

For multimode, GPU computations (with Nvidia cuda) is highly recommended. I have written a lot of cuda files to speed up simulations. It is controlled by "sim.gpu_yes=true or false".

For details, please see our paper: (link here; wait to be published).

There is a readme.pdf in the Documentation/ folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.
