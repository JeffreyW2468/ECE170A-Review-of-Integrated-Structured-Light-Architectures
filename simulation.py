import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to compute Stokes parameters
def compute_stokes(a1, a2, delta):
    s0 = a1**2 + a2**2
    s1 = a1**2 - a2**2
    s2 = 2 * a1 * a2 * np.cos(delta)
    s3 = 2 * a1 * a2 * np.sin(delta)
    return s0, s1, s2, s3

def compute_stokes_linear(a1, a2, delta):
    s0 = a1**2 + a2**2
    s1 = a1**2 - a2**2
    s2 = 2 * a1 * a2 * np.cos(delta)
    s3 = 0
    return s0, s1, s2, s3

# Normalize Stokes parameters
def normalize_stokes(s0, s1, s2, s3):
    return s1 / s0, s2 / s0, s3 / s0

# Generate points for Fig. 1a: Phase variation for channels 1, 3, 5
# Models alternating linear polarization states confined to the equatorial plane.
def generate_fig1a_points(num_points):
    points = []
    for i in range(num_points):
        a1 = 1
        a2 = 1
        delta = 2* np.pi * i / num_points  # Phase varies gradually from 0 to 2pi for channels 1, 3, 5
        points.append(compute_stokes_linear(a1, a2, delta))
    return points
  
# Generate points for Fig. 1b: Asymmetric linear polarization with amplitude imbalance
# Models elliptical polarization states due to asymmetry between horizontal and vertical components.
def generate_fig1b_points(num_points):
    points = []
    for i in range(num_points):
        a1 = 1
        a2 = 0.9 # Asymmetric y-component
        delta = 2 * np.pi * i / num_points + 0.1 * np.sin(2 * np.pi * i / num_points)
        points.append(compute_stokes(a1, a2, delta))
    return points

# Generate points for Fig. 1c: Asymmetric circular polarization
# Introduces sinusoidal amplitude variation to model transitions
# between linear, elliptical, and circular polarization states.
def generate_fig1c_points(num_points):
    points = []
    for i in range(num_points):
        a1 = 1
        a2 = 0.8 + 0.2 * np.sin(2 * np.pi * i / num_points)  # Slight amplitude variation
        delta = 2 * np.pi * i / num_points  # Phase variation
        points.append(compute_stokes(a1, a2, delta))
    return points

# Generate points for each case
num_points = 100
fig1a_points = generate_fig1a_points(num_points)
fig1b_points = generate_fig1b_points(num_points)
fig1c_points = generate_fig1c_points(num_points)

# Extract normalized Stokes parameters for plotting
S1_a, S2_a, S3_a = zip(*[normalize_stokes(s0, s1, s2, s3) for s0, s1, s2, s3 in fig1a_points])
S1_b, S2_b, S3_b = zip(*[normalize_stokes(s0, s1, s2, s3) for s0, s1, s2, s3 in fig1b_points])
S1_c, S2_c, S3_c = zip(*[normalize_stokes(s0, s1, s2, s3) for s0, s1, s2, s3 in fig1c_points])

# Plot the Poincaré sphere with calculated points
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the unit sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='lightblue', alpha=0.2, edgecolor='gray')

# Plot points for each case
ax.scatter(S1_a, S2_a, S3_a, color='red', s=20, label="Fig. 1a: Alternating Linear")
ax.scatter(S1_b, S2_b, S3_b, color='green', s=20, label="Fig. 1b: Asymmetric Linear")
ax.scatter(S1_c, S2_c, S3_c, color='blue', s=20, label="Fig. 1c: Asymmetric Circular")

ax.set_xlabel("$S_1$")
ax.set_ylabel("$S_2$")
ax.set_zlabel("$S_3$")
ax.set_title("Poincaré Sphere: Polarization Cases for Fig. 1")
ax.legend()
plt.show()
