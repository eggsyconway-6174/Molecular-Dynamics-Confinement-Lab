from vpython import *
import numpy as np
import math
import random
import time
import csv

# --- SHARED SIMULATION PARAMETERS ---
# N will be set by the loop
DIAMETER = 1.0          # Monomer diameter
MASS = 1.0              # Monomer mass
DT = 0.001              # Time step
N_STEPS = 10000         # Total simulation steps
EQUILIBRATION_STEPS = 2000 # Steps to ignore for averaging
TEMPERATURE = 10.0      # System temperature
K_HARMONIC = 100.0      # Bond spring constant
L_HARMONIC = DIAMETER + 0.1 # Bond equilibrium length
EPSILON = 1.0           # LJ potential depth
SIGMA = DIAMETER        # LJ potential size

# --- N values to test ---
N_values = range(20, 101, 10) # From 20 to 100, in steps of 10
N_BASE = 20.0 # Base N for scaling

# --- Derived Constants ---
SIGMA_SQ = SIGMA**2
LJ_CUTOFF = (2**(1.0/6.0)) * SIGMA
LJ_CUTOFF_SQ = LJ_CUTOFF**2

# --- CONFINEMENT PARAMETERS ---
R_SPHERE_BASE = 4.5  # Radius for N=20
EPSILON_WALL = 1.0
SIGMA_WALL = DIAMETER
WALL_CUTOFF = (2**(1.0/6.0)) * SIGMA_WALL
# R_SPHERE will be calculated inside the loop

# --- PBC PARAMETERS ---
# BOX_L will be calculated inside the loop

print(f"--- Starting Comparative Simulation ---")
print(f"Testing N values: {list(N_values)}")
print(f"Temperature T = {TEMPERATURE}")
print(f"Steps per N = {N_STEPS}")
print("---------------------------------------")


# --- VPYTHON SETUP (Create ONCE) ---
# Scene 1: Confined
scene_confined = canvas(title='Simulation 1: Spherical Confinement',
                        width=600, height=500, x=0, y=0,
                        center=vector(0,0,0), background=color.black)
# Create container, but make it invisible for now
container = sphere(canvas=scene_confined, pos=vector(0,0,0), radius=1, 
                   opacity=0.2, color=color.gray(0.5), visible=False)

# Scene 2: PBC
scene_pbc = canvas(title='Simulation 2: Periodic Boundary Conditions (PBC)',
                   width=600, height=500, x=600, y=0,
                   center=vector(0,0,0), background=color.black)
# Create box, but make it invisible for now
box_vis = box(canvas=scene_pbc, pos=vector(0,0,0),
              size=vector(1,1,1),
              color=color.gray(0.5), opacity=0.1, visible=False)

# --- Graphs (Shared, Create ONCE) ---
g_energies = graph(width=600, height=300, title='<b>System Energies</b>',
                   xtitle='<i>Time Step</i>', ytitle='<i>Energy</i>',
                   x=0, y=500, fast=False)
g_pe_c = gcurve(graph=g_energies, color=color.blue, label='PE (Confined)')
g_pe_p = gcurve(graph=g_energies, color=color.cyan, label='PE (PBC)')
g_ke_c = gcurve(graph=g_energies, color=color.red, label='KE (Confined)')
g_ke_p = gcurve(graph=g_energies, color=color.orange, label='KE (PBC)')

g_polymer = graph(width=600, height=300, title='<b>Polymer Properties</b>',
                  xtitle='<i>Time Step</i>', ytitle='<i>Distance</i>',
                  x=600, y=500, fast=False)
g_rg_c = gcurve(graph=g_polymer, color=color.green, label='Rg (Confined)')
g_rg_p = gcurve(graph=g_polymer, color=color.white, label='Rg (PBC)')
g_ree_c = gcurve(graph=g_polymer, color=color.purple, label='Ree (Confined)')
g_ree_p = gcurve(graph=g_polymer, color=color.yellow, label='Ree (PBC)')


# --- HELPER FUNCTIONS ---

def random_unit_vector():
    """Generates a random 3D unit vector."""
    phi = random.uniform(0, 2 * math.pi)
    costheta = random.uniform(-1, 1)
    theta = math.acos(costheta)
    x = math.sin(theta) * math.cos(phi)
    y = math.sin(theta) * math.sin(phi)
    z = math.cos(theta)
    return vector(x, y, z)

def minimum_image_vector(dr, L, half_L):
    """Applies the Minimum Image Convention (MIC) to a distance vector."""
    if dr.x > half_L:    dr.x -= L
    elif dr.x < -half_L: dr.x += L
    if dr.y > half_L:    dr.y -= L
    elif dr.y < -half_L: dr.y += L
    if dr.z > half_L:    dr.z -= L
    elif dr.z < -half_L: dr.z += L
    return dr

# --- PHYSICS: FORCE CALCULATIONS ---

def compute_forces_confined(particle_list, N, R_SPHERE):
    """Computes total forces for the confined system."""
    pe_lj = 0.0
    pe_h = 0.0
    pe_conf = 0.0
    
    for p in particle_list:
        p.force = vector(0,0,0)

    # --- Harmonic Bond Forces ---
    for i in range(N - 1):
        p_i = particle_list[i]
        p_j = particle_list[i+1]
        
        rij_vec = p_i.pos - p_j.pos
        r = rij_vec.mag
        
        if r < 1e-6: continue # Avoid division by zero
            
        force_mag_h = -K_HARMONIC * (r - L_HARMONIC)
        force_vec_h = (force_mag_h / r) * rij_vec
        
        p_i.force += force_vec_h
        p_j.force -= force_vec_h
        pe_h += 0.5 * K_HARMONIC * (r - L_HARMONIC)**2

    # --- Lennard-Jones Forces (Non-bonded) ---
    for i in range(N):
        for j in range(i + 2, N):  # Skip adjacent particles
            p_i = particle_list[i]
            p_j = particle_list[j]
            
            rij_vec = p_i.pos - p_j.pos
            r_sq = rij_vec.mag2
            
            if r_sq > LJ_CUTOFF_SQ:
                continue
                
            sr2 = SIGMA_SQ / r_sq
            sr6 = sr2**3
            sr12 = sr6**2
            
            force_mag_lj = 24 * EPSILON / r_sq * (2 * sr12 - sr6)
            force_vec_lj = force_mag_lj * rij_vec
            
            p_i.force += force_vec_lj
            p_j.force -= force_vec_lj
            pe_lj += 4 * EPSILON * (sr12 - sr6) # Un-shifted WCA
        
    # --- Spherical Confinement Force ---
    for p in particle_list:
        r_mag = p.pos.mag
        dist_to_wall = R_SPHERE - r_mag
        
        if dist_to_wall <= WALL_CUTOFF:
            r_sq = dist_to_wall**2
            if r_sq < 1e-10: r_sq = 1e-10 # Avoid explosion
            
            sr2 = (SIGMA_WALL**2) / r_sq
            sr6 = sr2**3
            sr12 = sr6**2
            
            force_mag = 24 * EPSILON_WALL / r_sq * (2 * sr12 - sr6)
            
            p.force -= force_mag * p.pos.norm() 
            pe_conf += 4 * EPSILON_WALL * (sr12 - sr6) # Un-shifted WCA
            
    return pe_lj, pe_h, pe_conf

def compute_forces_pbc(particle_list, N, L, half_L):
    """Computes total forces for the PBC system."""
    pe_lj = 0.0
    pe_h = 0.0
    
    for p in particle_list:
        p.force = vector(0,0,0)

    # --- Harmonic Bond Forces (with MIC) ---
    for i in range(N - 1):
        p_i = particle_list[i]
        p_j = particle_list[i+1]
        
        rij_vec = p_i.pos - p_j.pos
        rij_vec_mic = minimum_image_vector(rij_vec, L, half_L)
        r = rij_vec_mic.mag
        
        if r < 1e-6: continue
            
        force_mag_h = -K_HARMONIC * (r - L_HARMONIC)
        force_vec_h = (force_mag_h / r) * rij_vec_mic # Use MIC vector for force
        
        p_i.force += force_vec_h
        p_j.force -= force_vec_h
        pe_h += 0.5 * K_HARMONIC * (r - L_HARMONIC)**2

    # --- Lennard-Jones Forces (Non-bonded, with MIC) ---
    for i in range(N):
        for j in range(i + 2, N):
            p_i = particle_list[i]
            p_j = particle_list[j]
            
            rij_vec = p_i.pos - p_j.pos
            rij_vec_mic = minimum_image_vector(rij_vec, L, half_L)
            r_sq = rij_vec_mic.mag2
            
            if r_sq > LJ_CUTOFF_SQ:
                continue
                
            sr2 = SIGMA_SQ / r_sq
            sr6 = sr2**3
            sr12 = sr6**2
            
            force_mag_lj = 24 * EPSILON / r_sq * (2 * sr12 - sr6)
            force_vec_lj = force_mag_lj * rij_vec_mic # Use MIC vector
            
            p_i.force += force_vec_lj
            p_j.force -= force_vec_lj
            pe_lj += 4 * EPSILON * (sr12 - sr6) # Un-shifted WCA
            
    return pe_lj, pe_h

# --- PARTICLE INITIALIZATION ---

def initialize_particles(N, vis_scene, start_pos, max_dist_from_origin):
    """
    Creates a list of particle objects using a random walk.
    Ensures particles stay within bounds (sphere or box).
    """
    particles = []
    max_attempts = 1000
    
    for i in range(N):
        p_color = color.cyan
        if i == 0: p_color = color.red
        if i == N-1: p_color = color.yellow
        
        # Create particle in the correct scene
        p = sphere(canvas=vis_scene, radius=DIAMETER/2, color=p_color, make_trail=False)
        
        attempts = 0
        placed = False
        
        while attempts < max_attempts and not placed: 
            attempts += 1
            
            if i == 0:
                new_pos = start_pos
            else:
                new_pos = particles[i-1].pos + random_unit_vector() * L_HARMONIC
            
            # Check boundary
            if vis_scene == scene_confined:
                # Use magnitude for sphere
                if (new_pos.mag + DIAMETER/2) > max_dist_from_origin:
                    continue # Try a new random vector
            else:
                # Simple box check for initialization
                if abs(new_pos.x - start_pos.x) > max_dist_from_origin / 2 or \
                   abs(new_pos.y - start_pos.y) > max_dist_from_origin / 2 or \
                   abs(new_pos.z - start_pos.z) > max_dist_from_origin / 2:
                    continue

            # Check overlap
            valid = True
            for j in range(i):
                dr = new_pos - particles[j].pos
                # For PBC, this overlap check is not rigorous,
                # but it's okay for initialization.
                if dr.mag2 < DIAMETER**2:
                    valid = False
                    break
            
            if valid:
                p.pos = new_pos
                placed = True
                break
        
        if not placed:
            print(f"Warning: Could not place particle {i} without overlap. Placing anyway.")
            p.pos = new_pos
        
        p.mass = MASS
        p.velocity = vector(0, 0, 0)
        p.force = vector(0, 0, 0)
        particles.append(p)

    # Connect particles with visual bonds
    bonds = []
    for i in range(N - 1):
        bond = cylinder(canvas=vis_scene, pos=particles[i].pos, 
                        axis=particles[i+1].pos - particles[i].pos,
                        radius=0.05, color=color.white)
        bonds.append(bond)
        
    return particles, bonds

def initialize_velocities(particle_list, T, M):
    """Initializes velocities from Maxwell-Boltzmann distribution."""
    v_scale = np.sqrt(T / M)
    for p in particle_list:
        vx = random.gauss(0, v_scale)
        vy = random.gauss(0, v_scale)
        vz = random.gauss(0, v_scale)
        p.velocity = vector(vx, vy, vz)
    
    # Remove center-of-mass motion
    total_momentum = vector(0, 0, 0)
    for p in particle_list:
        total_momentum += p.mass * p.velocity
    avg_velocity = total_momentum / (len(particle_list) * M)
    for p in particle_list:
        p.velocity -= avg_velocity

# --- OBSERVABLE CALCULATIONS ---

def calculate_observables(particle_list, N, L, half_L, use_pbc=False):
    """Calculates KE, Rg, and Ree."""
    ke = 0.0
    cm = vector(0, 0, 0)
    
    # Need to "unwrap" coordinates for Rg and Ree in PBC
    unwrapped_pos = [vector(0,0,0)] * N
    unwrapped_pos[0] = particle_list[0].pos
    
    for p in particle_list:
        ke += 0.5 * p.mass * p.velocity.mag2
    
    if use_pbc:
        # Reconstruct "unwrapped" chain for Rg/Ree
        for i in range(1, N):
            dr = particle_list[i].pos - particle_list[i-1].pos
            # Reverse the MIC
            if dr.x > half_L:    dr.x -= L
            elif dr.x < -half_L: dr.x += L
            if dr.y > half_L:    dr.y -= L
            elif dr.y < -half_L: dr.y += L
            if dr.z > half_L:    dr.z -= L
            elif dr.z < -half_L: dr.z += L
            unwrapped_pos[i] = unwrapped_pos[i-1] + dr
        
        # Calculate CM from unwrapped positions
        for pos in unwrapped_pos:
            cm += pos
        cm /= N
        
        # Calculate Rg from unwrapped positions
        rg_sq = 0.0
        for pos in unwrapped_pos:
            rg_sq += (pos - cm).mag2
        
        if rg_sq < 0: rg_sq = 0 # Safety check
        rg = math.sqrt(rg_sq / N)
        ree = (unwrapped_pos[N-1] - unwrapped_pos[0]).mag
        
    else: # Confinement
        for p in particle_list:
            cm += p.pos
        cm /= N
        
        rg_sq = 0.0
        for p in particle_list:
            rg_sq += (p.pos - cm).mag2
        
        if rg_sq < 0: rg_sq = 0 # Safety check
        rg = math.sqrt(rg_sq / N)
        ree = (particle_list[N-1].pos - particle_list[0].pos).mag

    return ke, rg, ree

# --- NEW HELPER FUNCTION FOR CSV SAVING ---
def save_time_series_to_csv(filename, headers, time_col, data_cols):
    """Saves time series data to a CSV file."""
    print(f"Saving time series data to {filename}...")
    try:
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            # Combine time column with other data columns
            all_data = [time_col] + data_cols
            # Zip to transpose data from columns to rows
            rows = zip(*all_data)
            writer.writerows(rows)
        print(f"Successfully saved {filename}")
    except Exception as e:
        print(f"Error saving {filename}: {e}")

# --- MAIN EXPERIMENT LOOP ---

# Create a file to save the final averaged results
results_file = open("comparison_results_by_N.txt", "w")
results_file.write(f"N\tAvg_KE_C\tAvg_PE_C\tAvg_Rg_C\tAvg_Ree_C\tAvg_KE_P\tAvg_PE_P\tAvg_Rg_P\tAvg_Ree_P\n")
results_file.flush()


for N in N_values:
    print(f"\n=========================================")
    print(f"      Running Simulation for N = {N}")
    print(f"=========================================\n")

    # --- 1. Calculate N-dependent parameters ---
    R_SPHERE = R_SPHERE_BASE * (N / N_BASE)**(1.0/3.0)
    SPHERE_VOLUME = (4.0/3.0) * math.pi * R_SPHERE**3
    BOX_L = SPHERE_VOLUME**(1.0/3.0)
    BOX_HALF_L = BOX_L / 2.0
    
    print(f"--- System Parameters for N={N} ---")
    print(f"Sphere Radius R = {R_SPHERE:.3f}")
    print(f"Box Length L = {BOX_L:.3f} (Volume = {SPHERE_VOLUME:.3f})")
    print("---------------------------------")

    # Data lists for time series CSV
    time_data = []
    ke_c_data, pe_c_data, rg_c_data, ree_c_data = [], [], [], []
    ke_p_data, pe_p_data, rg_p_data, ree_p_data = [], [], [], []

    # --- 2. Update VPython Scenes ---
    g_energies.title = f'<b>System Energies (N={N})</b>'
    g_polymer.title = f'<b>Polymer Properties (N={N})</b>'
    # Clear data from previous run
    g_pe_c.data = []
    g_pe_p.data = []
    g_ke_c.data = []
    g_ke_p.data = []
    g_rg_c.data = []
    g_rg_p.data = []
    g_ree_c.data = []
    g_ree_p.data = []
    
    # Update confinement sphere
    container.radius = R_SPHERE
    container.visible = True
    scene_confined.range = R_SPHERE * 1.2
    scene_confined.center = vector(0,0,0)
    
    # Update PBC box
    box_vis.size = vector(BOX_L, BOX_L, BOX_L)
    box_vis.pos = vector(BOX_L/2, BOX_L/2, BOX_L/2)
    box_vis.visible = True
    scene_pbc.range = BOX_L * 0.7
    scene_pbc.center = vector(BOX_L/2, BOX_L/2, BOX_L/2)
    
    # --- 3. Initialize Systems for this N ---
    particles_c, bonds_c = initialize_particles(N, scene_confined, vector(0,0,0), R_SPHERE - DIAMETER/2)
    initialize_velocities(particles_c, TEMPERATURE, MASS)

    particles_p, bonds_p = initialize_particles(N, scene_pbc, vector(BOX_L/2, BOX_L/2, BOX_L/2), BOX_L)
    initialize_velocities(particles_p, TEMPERATURE, MASS)

    # --- 4. Calculate initial forces ---
    pe_lj_c, pe_h_c, pe_conf_c = compute_forces_confined(particles_c, N, R_SPHERE)
    pe_total_c = pe_lj_c + pe_h_c + pe_conf_c

    pe_lj_p, pe_h_p = compute_forces_pbc(particles_p, N, BOX_L, BOX_HALF_L)
    pe_total_p = pe_lj_p + pe_h_p

    # --- 5. Averages accumulators ---
    avg_ke_c, avg_pe_c, avg_rg_c, avg_ree_c = 0.0, 0.0, 0.0, 0.0
    avg_ke_p, avg_pe_p, avg_rg_p, avg_ree_p = 0.0, 0.0, 0.0, 0.0
    measurement_count = 0

    # --- 6. Main Simulation Loop for this N ---
    for step in range(N_STEPS):
        rate(500) # Keep visualization fast
        
        # --- A. CONFINED SIMULATION STEP ---
        old_forces_c = [p.force for p in particles_c]
        
        # Verlet 1: Update positions
        for p in particles_c:
            p.pos += p.velocity * DT + 0.5 * (p.force / p.mass) * DT**2
        
        # Update visual bonds
        for i in range(N - 1):
            bonds_c[i].pos = particles_c[i].pos
            bonds_c[i].axis = particles_c[i+1].pos - particles_c[i].pos

        # Calculate new forces
        pe_lj_c, pe_h_c, pe_conf_c = compute_forces_confined(particles_c, N, R_SPHERE)
        pe_total_c = pe_lj_c + pe_h_c + pe_conf_c
        
        # Verlet 2: Update velocities
        for i, p in enumerate(particles_c):
            p.velocity += 0.5 * ((old_forces_c[i] + p.force) / p.mass) * DT
            
        
        # --- B. PBC SIMULATION STEP ---
        old_forces_p = [p.force for p in particles_p]

        # Verlet 1: Update positions
        for p in particles_p:
            p.pos += p.velocity * DT + 0.5 * (p.force / p.mass) * DT**2
            
            # Apply PBC Wrap
            p.pos.x = p.pos.x % BOX_L
            p.pos.y = p.pos.y % BOX_L
            p.pos.z = p.pos.z % BOX_L
                
        # Update visual bonds (with MIC)
        for i in range(N - 1):
            p1 = particles_p[i].pos
            p2 = particles_p[i+1].pos
            dr_bond = p2 - p1
            dr_bond_mic = minimum_image_vector(dr_bond, BOX_L, BOX_HALF_L)
            bonds_p[i].pos = p1
            bonds_p[i].axis = dr_bond_mic

        # Calculate new forces
        pe_lj_p, pe_h_p = compute_forces_pbc(particles_p, N, BOX_L, BOX_HALF_L)
        pe_total_p = pe_lj_p + pe_h_p
        
        # Verlet 2: Update velocities
        for i, p in enumerate(particles_p):
            p.velocity += 0.5 * ((old_forces_p[i] + p.force) / p.mass) * DT
            
            
        # --- C. MEASURE & PLOT (every 10 steps) ---
        if step % 10 == 0:
            ke_c, rg_c, ree_c = calculate_observables(particles_c, N, BOX_L, BOX_HALF_L, use_pbc=False)
            ke_p, rg_p, ree_p = calculate_observables(particles_p, N, BOX_L, BOX_HALF_L, use_pbc=True)
            
            # Plot
            g_pe_c.plot(step, pe_total_c)
            g_pe_p.plot(step, pe_total_p)
            g_ke_c.plot(step, ke_c)
            g_ke_p.plot(step, ke_p)
            
            g_rg_c.plot(step, rg_c)
            g_rg_p.plot(step, rg_p)
            g_ree_c.plot(step, ree_c)
            g_ree_p.plot(step, ree_p)
            
            # Accumulate for averages (after equilibration)
            if step > EQUILIBRATION_STEPS:
                avg_ke_c += ke_c
                avg_pe_c += pe_total_c
                avg_rg_c += rg_c
                avg_ree_c += ree_c
                
                avg_ke_p += ke_p
                avg_pe_p += pe_total_p
                avg_rg_p += rg_p
                avg_ree_p += ree_p
                
                measurement_count += 1
            
            # Store data for time series
            time_data.append(step)
            ke_c_data.append(ke_c)
            pe_c_data.append(pe_total_c)
            rg_c_data.append(rg_c)
            ree_c_data.append(ree_c)
            ke_p_data.append(ke_p)
            pe_p_data.append(pe_total_p)
            rg_p_data.append(rg_p)
            ree_p_data.append(ree_p)
                
            # Check for explosion
            if not (np.isfinite(ke_c) and np.isfinite(ke_p)):
                print(f"ERROR: Simulation exploded at step {step}. Stopping.")
                break
    
    # --- 7. Save Time Series Data to CSV ---
    csv_headers_c = ['TimeStep', 'KE_Confined', 'PE_Confined', 'Rg_Confined', 'Ree_Confined']
    csv_data_c = [ke_c_data, pe_c_data, rg_c_data, ree_c_data]
    save_time_series_to_csv(f"timeseries_N_{N}_confined.csv", csv_headers_c, time_data, csv_data_c)
    
    csv_headers_p = ['TimeStep', 'KE_PBC', 'PE_PBC', 'Rg_PBC', 'Ree_PBC']
    csv_data_p = [ke_p_data, pe_p_data, rg_p_data, ree_p_data]
    save_time_series_to_csv(f"timeseries_N_{N}_pbc.csv", csv_headers_p, time_data, csv_data_p)

    # --- 8. Print Final Averages for this N ---
    print(f"\n--- Simulation Finished for N = {N} ---")
    if measurement_count > 0:
        avg_ke_c /= measurement_count
        avg_pe_c /= measurement_count
        avg_rg_c /= measurement_count
        avg_ree_c /= measurement_count
        
        avg_ke_p /= measurement_count
        avg_pe_p /= measurement_count
        avg_rg_p /= measurement_count
        avg_ree_p /= measurement_count

        print("\n--- Averages (post-equilibration) ---")
        print(f"| Observable | Confined (Sim 1) | PBC (Sim 2) |")
        print(f"|------------|--------------------|-------------|")
        print(f"| Avg. KE    | {avg_ke_c:18.4f} | {avg_ke_p:11.4f} |")
        print(f"| Avg. PE    | {avg_pe_c:18.4f} | {avg_pe_p:11.4f} |")
        print(f"| Avg. Rg    | {avg_rg_c:18.4f} | {avg_rg_p:11.4f} |")
        print(f"| Avg. Ree   | {avg_ree_c:18.4f} | {avg_ree_p:11.4f} |")
        
        # Write to results file
        results_file.write(f"{N}\t{avg_ke_c:.4f}\t{avg_pe_c:.4f}\t{avg_rg_c:.4f}\t{avg_ree_c:.4f}\t")
        results_file.write(f"{avg_ke_p:.4f}\t{avg_pe_p:.4f}\t{avg_rg_p:.4f}\t{avg_ree_p:.4f}\n")
        results_file.flush()
    else:
        print("No measurements taken post-equilibration.")

    print("---------------------------------")
    
    # --- 9. Clean up VPython objects before next loop ---
    container.visible = False
    box_vis.visible = False
    for p in particles_c: p.visible = False
    for b in bonds_c: b.visible = False
    for p in particles_p: p.visible = False
    for b in bonds_p: b.visible = False
    
    del particles_c, bonds_c, particles_p, bonds_p
    
    time.sleep(1) # Pause briefly

# --- End of N_values loop ---
results_file.close()
print("\n=========================================")
print("      All Simulations Finished")
print(f"Final results saved to 'comparison_results_by_N.txt'")
print("=========================================")
