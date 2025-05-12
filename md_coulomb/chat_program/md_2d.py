import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML

class MDSimulation2D:
    """
    A class for 2D molecular dynamics simulations with a radial inter-particle potential.
    Supports NVE (microcanonical) and NVT (canonical) ensembles, with periodic or reflecting boundaries.
    The default potential is repulsive Coulomb (V(r) = 1/r), but a custom potential function can be provided.
    """
    def __init__(self, N, L, dt, steps, 
                 initial_temperature=1.0, target_temperature=None,
                 ensemble='NVE', thermostat='rescale',
                 potential_func=None, 
                 boundary='periodic', init_positions='random', seed=None):
        """
        Initialize the simulation.
        Parameters:
            N (int): Number of particles.
            L (float or tuple): Box size. If float, uses square box [0, L] in both x and y.
            dt (float): Time step for integration.
            steps (int): Number of simulation steps to run.
            initial_temperature (float): Initial temperature for velocity initialization.
            target_temperature (float or None): Target temperature for NVT ensemble (defaults to initial_temperature).
            ensemble (str): 'NVE' or 'NVT' ensemble.
            thermostat (str): Thermostat method for NVT ('rescale' or 'langevin'). Only used if ensemble='NVT'.
            potential_func (callable or None): Custom radial potential function V(r). 
                If None, uses default repulsive Coulomb potential V(r) = 1/r.
            boundary (str): 'periodic' or 'reflect' boundary conditions.
            init_positions (str): 'random' or 'grid' initial placement of particles.
            seed (int or None): Random seed for reproducibility.
        """
        # Simulation parameters
        self.N = N
        # Box dimensions (Lx, Ly)
        if np.isscalar(L):
            self.Lx = self.Ly = float(L)
        else:
            self.Lx, self.Ly = float(L[0]), float(L[1])
        self.dt = dt
        self.steps = steps
        self.ensemble = ensemble
        self.thermostat = thermostat
        self.initial_temperature = initial_temperature
        # If target_temperature is not provided for NVT, use initial_temperature
        self.target_temperature = target_temperature if target_temperature is not None else initial_temperature
        self.boundary = boundary
        # Potential function: default 1/r Coulomb if not provided
        if potential_func is None:
            self.potential_func = lambda r: 1.0 / r  # Coulomb potential (k=1 for simplicity)
        else:
            self.potential_func = potential_func
        # Set random seed if provided
        if seed is not None:
            np.random.seed(seed)
        # Initialize particle positions and velocities
        self.mode = init_positions
        self.positions = self._initialize_positions(mode=init_positions)
        self.velocities = self._initialize_velocities(self.initial_temperature)
        # Lists to record observables over time
        self.kinetic_energy_history = []
        self.potential_energy_history = []
        self.total_energy_history = []
        self.temperature_history = []
        self.pressure_history = []

    def reset(self):
        self._initialize_positions(mode = self.mode)
        self._initialize_velocities(self.initial_temperature)
    
    def density(self):
        return self.N / (self.Lx * self.Ly)

    def _initialize_positions(self, mode='random'):
        """Initialize particle positions either on a square grid or randomly within the box."""
        if mode == 'grid':
            # Place particles on a grid inside the box
            n_side = int(np.ceil(np.sqrt(self.N)))
            xs = np.linspace(0, self.Lx, n_side, endpoint=False) + (self.Lx / n_side) / 2
            ys = np.linspace(0, self.Ly, n_side, endpoint=False) + (self.Ly / n_side) / 2
            xv, yv = np.meshgrid(xs, ys)
            grid_points = np.vstack([xv.ravel(), yv.ravel()]).T
            positions = grid_points[:self.N].copy()
        else:  # 'random' placement
            # Uniform random distribution in box
            positions = np.zeros((self.N, 2))
            positions[:, 0] = np.random.rand(self.N) * self.Lx
            positions[:, 1] = np.random.rand(self.N) * self.Ly
            # Optionally, avoid extremely close initial pairs by iterative adjustment
            for _ in range(5):  # attempt a few times
                too_close = False
                for i in range(self.N):
                    # Distance from particle i to all others
                    diff = positions[i] - positions  # differences to others
                    if self.boundary == 'periodic':
                        # Apply periodic wrapping for distance calculation
                        diff[:, 0] -= self.Lx * np.round(diff[:, 0] / self.Lx)
                        diff[:, 1] -= self.Ly * np.round(diff[:, 1] / self.Ly)
                    dist = np.sqrt(diff[:, 0]**2 + diff[:, 1]**2)
                    dist[i] = np.inf  # ignore self-distance
                    # If any particle is too close (e.g., closer than 10% of box min dimension), reassign its position
                    if np.any(dist < 0.1 * min(self.Lx, self.Ly)):
                        positions[i, 0] = np.random.rand() * self.Lx
                        positions[i, 1] = np.random.rand() * self.Ly
                        too_close = True
                if not too_close:
                    break
        return positions

    def _initialize_velocities(self, temperature):
        """
        Initialize velocities from a Maxwell-Boltzmann distribution (Gaussian) for a given temperature.
        Removes any net center-of-mass drift so total momentum is zero.
        """
        # For mass=1 and k_B=1, variance = temperature, so std dev = sqrt(T)
        sigma = np.sqrt(temperature)
        velocities = np.random.normal(0.0, sigma, size=(self.N, 2))
        # Remove center-of-mass velocity (net momentum)
        com_velocity = np.mean(velocities, axis=0)
        velocities -= com_velocity
        return velocities

    def _compute_forces_and_potential(self):
        """
        Compute inter-particle forces and total potential energy using the radial potential function.
        Uses the minimum image convention for periodic boundaries.
        Returns:
            forces (ndarray): Array of shape (N,2) of forces on each particle.
            potential_energy (float): Total potential energy of the system.
        """
        # Pairwise displacement vectors: diff[i,j] = r_i - r_j
        diff = self.positions[:, np.newaxis, :] - self.positions[np.newaxis, :, :]
        if self.boundary == 'periodic':
            # Apply minimum image convention
            diff[..., 0] -= self.Lx * np.round(diff[..., 0] / self.Lx)
            diff[..., 1] -= self.Ly * np.round(diff[..., 1] / self.Ly)
        # Compute distance magnitudes
        dist_sq = diff[..., 0]**2 + diff[..., 1]**2
        np.fill_diagonal(dist_sq, np.inf)  # prevent self-interaction (set self-distances to infinity)
        distances = np.sqrt(dist_sq)
        # Compute potential matrix V(r) for all pairs (i,j)
        # The potential function is assumed to accept numpy arrays.
        V_matrix = self.potential_func(distances)
        # Total potential energy = 1/2 sum_{i != j} V(r_ij) (each pair counted twice in matrix)
        potential_energy = 0.5 * np.nansum(V_matrix)
        # Compute force magnitudes via numerical differentiation: F(r) = -dV/dr.
        epsilon = 1e-8  # small increment for derivative
        finite_mask = np.isfinite(distances)  # mask for real distances (exclude infinities)
        F_mag = np.zeros_like(distances)
        if np.any(finite_mask):
            r_vals = distances[finite_mask]
            # Avoid negative or zero radii for differentiation
            r_plus = r_vals + epsilon
            r_minus = np.where(r_vals > epsilon, r_vals - epsilon, 1e-12)
            # Evaluate potential at r+ and r- for those distances
            V_plus = self.potential_func(r_plus)
            V_minus = self.potential_func(r_minus)
            # Derivative dV/dr ≈ (V_plus - V_minus) / (r_plus - r_minus)
            dV_dr = (V_plus - V_minus) / (r_plus - r_minus)
            F_mag[finite_mask] = -dV_dr  # negative derivative gives force magnitude
        # Calculate force vectors from magnitudes
        # Unit vector from j to i = diff_ij / r_ij
        unit_vec = np.zeros_like(diff)
        unit_vec[..., 0] = np.divide(diff[..., 0], distances, out=np.zeros_like(diff[..., 0]), where=np.isfinite(distances))
        unit_vec[..., 1] = np.divide(diff[..., 1], distances, out=np.zeros_like(diff[..., 1]), where=np.isfinite(distances))
        # Force on particle i due to j: F_ij = F_mag[i,j] * unit_vector[i,j]
        # Total force on particle i = sum_j F_ij
        forces = np.zeros((self.N, 2))
        forces[:, 0] = np.nansum(F_mag * unit_vec[..., 0], axis=1)
        forces[:, 1] = np.nansum(F_mag * unit_vec[..., 1], axis=1)
        return forces, potential_energy

    def _apply_boundary(self):
        """Apply boundary conditions to particle positions (and adjust velocities if reflecting)."""
        if self.boundary == 'periodic':
            # Wrap positions that leave the box (torus topology)
            self.positions[:, 0] = np.mod(self.positions[:, 0], self.Lx)
            self.positions[:, 1] = np.mod(self.positions[:, 1], self.Ly)
        elif self.boundary == 'reflect':
            # Reflective walls: if a particle crosses, reflect its position and invert velocity
            # X direction
            over_left = self.positions[:, 0] < 0
            over_right = self.positions[:, 0] > self.Lx
            self.positions[over_left, 0] = -self.positions[over_left, 0]
            self.positions[over_right, 0] = 2*self.Lx - self.positions[over_right, 0]
            self.velocities[over_left | over_right, 0] *= -1  # reverse velocity for reflected particles
            # Y direction
            over_bottom = self.positions[:, 1] < 0
            over_top = self.positions[:, 1] > self.Ly
            self.positions[over_bottom, 1] = -self.positions[over_bottom, 1]
            self.positions[over_top, 1] = 2*self.Ly - self.positions[over_top, 1]
            self.velocities[over_bottom | over_top, 1] *= -1  # reverse velocity for reflected particles

    def step(self):
        """
        Perform one integration time step using Velocity Verlet algorithm.
        Applies forces, updates positions and velocities, and handles thermostat if NVT.
        """
        # Compute forces and potential energy at current positions
        forces, potential_energy = self._compute_forces_and_potential()
        # Velocity Verlet integration:
        # 1. Half-step velocity update: v(t+0.5dt) = v(t) + 0.5 * a(t) * dt (with a = F/m, m=1)
        self.velocities += 0.5 * forces * self.dt
        # 2. Position update: r(t+dt) = r(t) + v(t+0.5dt) * dt
        self.positions += self.velocities * self.dt
        # 3. Apply boundary conditions (this may adjust positions and velocities if reflecting)
        self._apply_boundary()
        # 4. Compute forces at new positions (t+dt)
        new_forces, new_potential_energy = self._compute_forces_and_potential()
        # 5. Complete velocity update: v(t+dt) = v(t+0.5dt) + 0.5 * a(t+dt) * dt
        self.velocities += 0.5 * new_forces * self.dt
        # Use updated potential energy for this time step
        potential_energy = new_potential_energy
        # Compute kinetic energy (m=1 for all particles)
        kinetic_energy = 0.5 * np.sum(self.velocities**2)
        # Instantaneous temperature (using k_B=1): T = 2K / (dof).
        # In 2D with N particles, dof = 2N, so T = K / (N*0.5)? Actually T = (2K) / (2N) = K/N.
        temperature = kinetic_energy / self.N
        # Compute pressure via virial theorem (for 2D): P = [sum_i (m v_i^2) + sum_i (r_i · F_i)] / (2V)
        kinetic_term = np.sum(self.velocities**2)  # = 2 * kinetic_energy for m=1
        virial_term = np.sum(self.positions[:, 0] * new_forces[:, 0] + self.positions[:, 1] * new_forces[:, 1])
        pressure = (kinetic_term + virial_term) / (2 * self.Lx * self.Ly)
        # Thermostat for NVT ensemble
        if self.ensemble.upper() == 'NVT':
            if self.thermostat == 'rescale':
                # Velocity rescaling: scale velocities to match target temperature
                lambda_factor = np.sqrt(self.target_temperature / (temperature + 1e-16))
                self.velocities *= lambda_factor
                # Recalculate energies after rescaling to record the adjusted values
                kinetic_energy = 0.5 * np.sum(self.velocities**2)
                temperature = kinetic_energy / self.N
            elif self.thermostat == 'langevin':
                # Langevin thermostat: include a drag (gamma) and random force to mimic heat bath
                gamma = 0.1  # friction coefficient
                sigma_v = np.sqrt(2 * gamma * self.target_temperature * self.dt)
                # Apply friction (damp velocities) and random thermal kicks
                self.velocities *= (1 - gamma * self.dt)
                self.velocities += sigma_v * np.random.normal(0, 1, size=self.velocities.shape)
                # Recalculate kinetic energy and temperature after Langevin thermostat
                kinetic_energy = 0.5 * np.sum(self.velocities**2)
                temperature = kinetic_energy / self.N
        # Record observables
        self.kinetic_energy_history.append(kinetic_energy)
        self.potential_energy_history.append(potential_energy)
        self.total_energy_history.append(kinetic_energy + potential_energy)
        self.temperature_history.append(temperature)
        self.pressure_history.append(pressure)

    def run(self):
        """Run the simulation for the specified number of steps."""
        for _ in range(self.steps):
            self.step()

    def plot_energies(self):
        """Plot kinetic, potential, and total energy as a function of time step."""
        t = np.arange(len(self.total_energy_history))
        plt.figure(figsize=(6,4))
        plt.plot(t, self.kinetic_energy_history, label='Kinetic Energy')
        plt.plot(t, self.potential_energy_history, label='Potential Energy')
        plt.plot(t, self.total_energy_history, label='Total Energy', linestyle='--', linewidth=2)
        plt.xlabel('Step')
        plt.ylabel('Energy')
        plt.title('Energy vs. Time')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_temperature(self):
        """Plot temperature vs. time step, with a line for target temperature if NVT."""
        t = np.arange(len(self.temperature_history))
        plt.figure(figsize=(6,4))
        plt.plot(t, self.temperature_history, label='Temperature', color='orange')
        if self.ensemble.upper() == 'NVT':
            plt.axhline(self.target_temperature, color='gray', linestyle='--', label='Target T')
        plt.xlabel('Step')
        plt.ylabel('Temperature')
        plt.title('Temperature vs. Time')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_pressure(self):
        """Plot pressure vs. time step, with an ideal gas reference line (NkT/V)."""
        t = np.arange(len(self.pressure_history))
        plt.figure(figsize=(6,4))
        plt.plot(t, self.pressure_history, label='Pressure', color='green')
        # Ideal gas pressure for reference (NkT/V)
        ideal_P = (self.N * self.target_temperature) / (self.Lx * self.Ly)
        plt.axhline(ideal_P, color='gray', linestyle='--', label='Ideal Gas P')
        plt.xlabel('Step')
        plt.ylabel('Pressure')
        plt.title('Pressure vs. Time')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_pair_correlation(self, bins=50):
        """Compute and plot the radial pair correlation function g(r) from current positions."""
        # Compute all pair distances using current positions
        diff = self.positions[:, np.newaxis, :] - self.positions[np.newaxis, :, :]
        if self.boundary == 'periodic':
            # Apply periodic wrap for distances
            diff[..., 0] -= self.Lx * np.round(diff[..., 0] / self.Lx)
            diff[..., 1] -= self.Ly * np.round(diff[..., 1] / self.Ly)
        dist_matrix = np.sqrt(diff[..., 0]**2 + diff[..., 1]**2)
        # Exclude self-distances
        mask = ~np.eye(self.N, dtype=bool)
        pair_distances = dist_matrix[mask]
        # Histogram of pair distances (up to half the smallest box dimension)
        r_max = 0.5 * min(self.Lx, self.Ly)
        hist, edges = np.histogram(pair_distances, bins=bins, range=(0, r_max))
        # Compute radial distribution function g(r)
        area = self.Lx * self.Ly
        density = self.N / area
        dr = edges[1] - edges[0]
        r_centers = edges[:-1] + dr/2
        # Area of each annular shell (2D): π*(r_outer^2 - r_inner^2)
        shell_areas = np.pi * (edges[1:]**2 - edges[:-1]**2)
        # Each pair was counted twice in `pair_distances` (i->j and j->i), 
        # divide by N to get average per particle in shell, and normalize by area and density
        g_r = hist / (self.N * density * shell_areas)
        # Plot g(r)
        plt.figure(figsize=(6,4))
        plt.plot(r_centers, g_r, drawstyle='steps-mid')
        plt.xlabel('Distance r')
        plt.ylabel('g(r)')
        plt.title('Radial Pair Correlation Function')
        plt.tight_layout()
        plt.show()

    def plot_positions(self, show_vectors=False, vector_scale=0.1):
        """
        Plot particle positions in the 2D box. Optionally overlay velocity vectors.
        Parameters:
            show_vectors (bool): If True, draw arrows for velocity vectors.
            vector_scale (float): Scale factor for velocity arrows (smaller value = longer arrows).
        """
        plt.figure(figsize=(6,6))
        plt.scatter(self.positions[:, 0], self.positions[:, 1], s=50, color='blue', label='Particles')
        if show_vectors:
            # Draw velocity vectors as arrows
            plt.quiver(self.positions[:, 0], self.positions[:, 1],
                       self.velocities[:, 0], self.velocities[:, 1],
                       angles='xy', scale_units='xy', scale=1/vector_scale, color='red', width=0.003,
                       label='Velocity')
        plt.xlim(0, self.Lx)
        plt.ylim(0, self.Ly)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Particle Positions' + (' with Velocity Vectors' if show_vectors else ''))
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def animate_trajectories(self, interval=50, steps_per_frame=1, show_vectors=False, vector_scale=0.1, frames=None, save_path=None):

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly)
        ax.set_aspect('equal')
        ax.set_title('MD Simulation Animation')

        scatter = ax.scatter(self.positions[:, 0], self.positions[:, 1], s=50, color='blue')
        if show_vectors:
            quiver = ax.quiver(self.positions[:, 0], self.positions[:, 1],
                            np.zeros_like(self.velocities[:, 0]),
                            np.zeros_like(self.velocities[:, 1]),
                            angles='xy', scale_units='xy', scale=1/vector_scale, color='red', width=0.003)
        else:
            quiver = None

        def update(frame_num):
            if frame_num > 0:
                for _ in range(steps_per_frame):
                    self.step()
                scatter.set_offsets(self.positions)
            if quiver:
                quiver.set_offsets(self.positions)
                quiver.set_UVC(self.velocities[:, 0], self.velocities[:, 1])
                return scatter, quiver
            else:
                return (scatter,)
            

        if frames is None:
            frames = self.steps // steps_per_frame

        ani = animation.FuncAnimation(
            fig, update, frames=frames,
            interval=interval, blit=True
        )
        
        if save_path:
            if save_path.lower().endswith('.mp4'):
                pass
                writer = animation.FFMpegWriter(fps=1000//interval)
            elif save_path.lower().endswith('.gif'):
                pass
                writer = animation.PillowWriter(fps=1000//interval)
            else:
                raise ValueError("Unsupported file format. Use .mp4 or .gif")
            ani.save(save_path, writer=writer)
        else:
            return HTML(ani.to_jshtml())

