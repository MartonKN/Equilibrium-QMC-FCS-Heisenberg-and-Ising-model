# Equilibrium-QMC-FCS-Heisenberg-and-Ising-model

This folder contains the code we needed to determine the full counting statistics of the staggered magnetization in the Heisenberg model, https://en.wikipedia.org/wiki/Heisenberg_model_(quantum). 

Corresponding paper: "Experimental realization of a long-range antiferromagnet in the Hubbard model with ultracold atoms", arXiv:1612.08436

About the code:

We use a stochastic series expansion quantum Monte Carlo to determine thermal averages. We sample the spin states of the system by updating them using Markov chain Monte Carlo updates. These updates affect a large amount of spins at once, hence they are much more effective than simple single spin updates. More about this method here: http://physics.bu.edu/~sandvik/vietri/sse.pdf

The code can be compiled using the makefile. The way to run the code:

./equilibrium_QMC TEMPERATURE NUMBER_OF_INDEPENDENT_RUNS MONTE_CARLO_UPDATES_PER_RUN MONTE_CARLO_BURN_IN_STEPS_PER_RUN  > myfile.txt

TEMPERATURE = temperature of the spin model

NUMBER_OF_INDEPENDENT_RUNS = how many times do we restart the simulation from a random initial condition

MONTE_CARLO_UPDATES_PER_RUN = how many Monte Carlo updates (cluster spin updates) do we perform per run

MONTE_CARLO_BURN_IN_STEPS_PER_RUN = number of Monte Carlo steps at the beginning of the simulation, until we do not measure any physical quantity (the so-called burn-in period)


