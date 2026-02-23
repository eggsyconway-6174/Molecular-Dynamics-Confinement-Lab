# 🧪 The Polymer Playground: Confinement vs. Bulk Dynamics

Welcome to a 3D computational experiment where we force digital molecules to fight against the laws of entropy. This isn't just a script; it's a dual-universe simulation comparing how a polymer chain survives in a **Caged Sphere** versus an **Infinite Void**.

## 🚀 The Mission: Why Build This?
In biology and chemistry, polymers (like DNA or proteins) don't always live in wide-open spaces. They get squeezed into cell nuclei or pushed through tiny pores. I built this to visualize the "Squeeze Factor." 

What happens to the shape of a molecule when you steal its space? Does it coil like a spring or collapse like a string? We’re using **Molecular Dynamics (MD)** to find out.

---

## 🏗️ Inside the Engine: What’s Happening?

### 1. The "Laws of Nature" (The Force Fields)
To make the particles move realistically, the code calculates three types of "Push and Pull":
* **The Invisible Bones (Harmonic Bonds):** Every monomer is connected by a "Hookean Spring." If they get too far, they snap back; too close, they push away. 
* **The "Get Away From Me" Force (WCA Potential):** We use a version of the Lennard-Jones potential to ensure no two atoms ever overlap. They have "hard shells," just like real matter.
* **The Force Field (Spherical Confinement):** In Scene 1, there is a mathematical wall. If a particle touches it, it gets blasted back toward the center using an inverse-square repulsion.



### 2. The Time Traveler (Velocity Verlet Integrator)
Standard "speed = distance/time" isn't enough for physics. The code uses the **Velocity Verlet Algorithm**. It predicts where a particle will be, calculates the force at that new spot, and then corrects the velocity. It’s a self-correcting loop that keeps the energy stable so the simulation doesn't "explode" into infinity.

### 3. The Infinite Loop (Minimum Image Convention)
In Scene 2, we use **PBC (Periodic Boundary Conditions)**. When a particle exits the right side of the box, it instantly "teleports" back in from the left. To make the physics work, the code uses **MIC (Minimum Image Convention)**—particles only "see" the version of their neighbors that is closest to them, even if that neighbor is technically in a ghost-universe next door.



---

## 📊 The "Aha!" Moments (Outcomes)

As the simulation runs, it generates a `results_by_N.txt` file. This is where the real science happens. We track:

* **Radius of Gyration ($R_g$):** The "fluffiness" of the polymer. In the sphere, you’ll notice $R_g$ plateaus because the walls won't let it grow. In the box, it expands freely.
* **End-to-End Distance ($R_{ee}$):** The distance from the "head" to the "tail." It’s a measure of how stretched the chain is.
* **The Scaling Law:** By running $N=20$ to $N=100$, we can actually calculate the **Flory Exponent**, proving that the polymer's size scales with its length in a very specific, predictable way.



---

## 🧰 The Tech Stack
* **VPython:** For the 3D "Neon-Cyber" real-time rendering.
* **NumPy:** To handle the heavy-lifting of the vector mathematics.
* **CSV/File IO:** To export our raw data for the final research report.

## 🧬 The "Fun" Part
The coolest bit of code is the **Coordinate Unwrapper**. When the polymer in the PBC box wraps around the edge, it looks broken to the eye. But the code "unwraps" it in the background, treating it like a continuous silk thread so the math for $R_g$ stays perfect. It’s like seeing the code and the reality at the same time.

---
*Built to explore the statistical mechanics of soft matter... and because watching 3D atoms bounce around is incredibly satisfying.*
