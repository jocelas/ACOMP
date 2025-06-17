def MC_step(x, beta, displacement):
    """Perform a single Monte Carlo step based of precalculated displacement

    Args:
        x (float): Current position
        T (float): Temperature
        displacement (float): Displacement to apply

    Returns:
        tuple: New position and whether the move was accepted    
    """
    x_new = x + displacement
    U_x = U(x)
    U_x_new = U(x_new)

    if U_x_new < U_x:
        return x_new, True
    else:
        p_accept = np.exp(-(U_x_new - U_x) * beta)
        if np.random.rand() < p_accept:
            return x_new, True
        else:
            return x, False
        

        
def Parallel_Tempering_step(x1, x2, beta1, beta2):
    """Perform a parallel tempering step between two temperatures by swapping positions if accepted

    Args:
        x1 (float): Current position at temperature T1
        x2 (float): Current position at temperature T2
        beta1 (float): Inverse temperature 1 (1/T1)
        beta2 (float): Inverse temperature 2 (1/T2)

    Returns:
        whether the swap was accepted
    """
    U_x1 = U(x1)
    U_x2 = U(x2)

    

    p_swap = np.exp((beta1 - beta2) * (U_x1 - U_x2))
    if np.random.rand() < p_swap:
        return x2, x1, True
    else:
        return x1, x2, False


def precalculate_displacements(maximum_displacement, n_steps):
    """Precalculate displacements for the Monte Carlo steps

    Args:
        maximum_displacement (float): Maximum displacement allowed
        n_steps (int): Number of steps to precalculate

    Returns:
        np.ndarray: Array of displacements
    """
    return np.random.rand(n_steps) * maximum_displacement * 2 - maximum_displacement

def simulation(n_steps, betas = betas, x_initial = x_initial, parallel_tempering_after = parallel_tempering_after, maximum_displacement = maximum_displacement, parallel_tempering = True):
    """Run the Monte Carlo simulation"""

    n_temperatures = len(betas)
    displacements = [precalculate_displacements(maximum_displacement, n_steps) for i in range(n_temperatures)]

    x = np.full(n_temperatures, x_initial)
    energies = np.zeros((n_temperatures, n_steps))
    mc_acceptance_ratio = np.zeros(n_temperatures)
    swap_acceptance_ratio = np.zeros((n_temperatures, n_temperatures))
    
    single_trajectory = np.empty(n_steps)

    for step in range(n_steps):
        for i in range(n_temperatures):
            x_new, accepted = MC_step(x[i], betas[i], displacements[i][step])
            energies[i, step] = U(x_new)
            if accepted:
                x[i] = x_new
                mc_acceptance_ratio[i] += 1

        if step % parallel_tempering_after == 0 and step > 0 and parallel_tempering:
            # Perform parallel tempering step
            for i in range(n_temperatures - 1):
                x[i], x[i + 1], swap_accepted = Parallel_Tempering_step(x[i], x[i + 1], betas[i], betas[i + 1])
                swap_acceptance_ratio[i, i + 1] += swap_accepted
                swap_acceptance_ratio[i + 1, i] += swap_accepted

        single_trajectory[step] = x[2]  # Store the trajectory of the first temperature

    mc_acceptance_ratio /= n_steps
    swap_acceptance_ratio /= (n_steps // parallel_tempering_after)

    return energies, mc_acceptance_ratio, swap_acceptance_ratio, single_trajectory


def trajectory_histogram(trajectory, n_bins=100):
    """Plot histogram of the trajectory

    Args:
        trajectory (np.ndarray): Trajectory data
        n_bins (int): Number of bins for the histogram
    """
    plt.figure(figsize=(8, 4))
    plt.hist(trajectory, bins=n_bins, density=True, alpha=0.6, color='g', label='Trajectory Histogram')

    # Plot the potential energy for reference
    x = np.linspace(-2.5, 2.5, 1000)
    U_x = np.vectorize(U)(x)
    plt.plot(x, U_x, 'r-', lw=2, label='Potential Energy')
    plt.title('Trajectory Histogram')
    plt.xlabel('Position')
    plt.ylabel('Density')
    plt.grid()
    plt.legend()
    plt.show()

