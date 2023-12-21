from cliffalg import Algebra

from math import pi
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


R3 = Algebra('space')
e_x, e_y, e_z, b_xy, b_xz, b_yz, t_xyz = R3.blades()[1:]

# Bivectors specifying the plane of rotation
rot_planes = [
    b_xy,
    b_xy+b_xz+b_yz,
    (3*b_xy + 4*b_xz)/5,
    
]   


wait = 25             # 25 ms -> 40 fps 
frames_per_plane = 400
frames = len(rot_planes)*(frames_per_plane)

# Make all 8 corners of a cube
vectors = [a*e_x + b*e_y + c*e_z for a in (-1, 1) for b in (-1, 1) for c in (-1,1)]

# A tracing of the corners that gets all edges.
trace = (0, 1, 3, 7, 5, 4, 6, 2, 0, 4, 5, 1, 3, 2, 6, 7)
arr = np.array([vectors[j].comps(1) for j in trace])


#This is a rotor sandwich product.
def R3_rotate(vec, plane, angle):
    return (-plane*angle/2).exp()*vec*(plane*angle/2).exp()

# Prepare plot
fig = plt.figure()
ax = plt.axes(projection='3d', proj_type='ortho')
ax.set_axis_off()
points, = ax.plot3D(arr[:, 0], arr[:, 1], arr[:, 2], 'o-')

# Plot plane of rotation
xx, yy = np.meshgrid((-0.3, 0.3), (-0.3, 0.3))
normal = (-t_xyz*rot_planes[0]).comps(1)           # Get components of normal vector
z = (-normal[0]*xx - normal[1]*yy)/normal[2]
plane = [ax.plot_surface(xx, yy, z, alpha=0.4)]

# Plot axis of rotation
axis, = ax.plot3D(np.array([0, normal[0]]), np.array([0, normal[1]]), np.array([0, normal[2]]), '-')

# Animation function
def animate(i, vectors, points, plane):
    global p
    if i % frames_per_plane == 0:    # Switch plane
        p = i//frames_per_plane
        normal = (-t_xyz*rot_planes[p]).comps(1)
        z = (-normal[0]*xx - normal[1]*yy)/normal[2]
        plane[0].remove()           # Delete old plane
        plane[0] = ax.plot_surface(xx, yy, z, alpha=0.4)

        axis.set_data(np.array([0, normal[0]]), np.array([0, normal[1]]))
        axis.set_3d_properties(np.array([0, normal[2]]), 'z')
        
    else:
        for i, old_vec in enumerate(vectors):
            # Rotate each vector with the rotor sandwich.
            vectors[i] = R3_rotate(old_vec, rot_planes[p], 2*pi/(frames_per_plane-1))
        
        arr = np.array([vectors[j].comps(1) for j in trace])    # Remake array

        points.set_data(arr[:, 0], arr[:, 1])
        points.set_3d_properties(arr[:, 2], 'z')


anim = FuncAnimation(fig, animate, frames=frames, fargs=(vectors, points, plane), interval=wait, repeat=True)
#anim.save('rotor.gif', fps=40)
plt.show()


