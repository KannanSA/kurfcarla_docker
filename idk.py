# OVITO Pro Python script for C60 nanocluster visualization
# -----------------------------------------------------------
#
# *** IMPORTANT ***
# This is a complete, standalone script.
#
# *** HOW TO USE ***
# 1. Make sure your 'C60.xyz' file is accessible at the path specified below.
# 2. In OVITO Pro, run this script from the top menu: File -> Run Python Script...
# 3. The script will automatically perform all steps: load the data, color the
#    bonds, and save a high-quality image and a rotation movie.

import ovito
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, PythonScriptModifier
from ovito.vis import TachyonRenderer, Viewport
from ovito.data import DataCollection
from ovito.pipeline import Pipeline, StaticSource
import numpy as np
import os

print("Running standalone C60 visualization script...")

# --- Python Script Function to Identify Bond Types and Set Colors ---
# This function is used by the PythonScriptModifier to perform the data analysis.
def setup_coloring(frame: int, data: DataCollection):
    """
    This function finds rings in the structure and directly creates 'Color'
    properties for both particles and bonds.
    """
    if not data.particles or not data.particles.bonds:
        return

    # Set a uniform color for all particles by creating a 'Color' property.
    num_particles = data.particles.count
    particle_color = (0.8, 0.8, 0.8) # Neutral grey
    data.particles_.create_property('Color', data=[particle_color] * num_particles)

    # Access the bond topology to perform ring-finding.
    bond_topology = data.particles.bonds.topology
    if bond_topology is None: return

    # Build a neighbor list for efficient traversal.
    neighbor_list = [[] for _ in range(num_particles)]
    for p1, p2 in bond_topology:
        neighbor_list[p1].append(p2)
        neighbor_list[p2].append(p1)

    # Find all bonds that are part of a 5-membered ring.
    pentagon_bonds = set()
    for i in range(num_particles):
        stack = [(i, [i])]
        while stack:
            current_atom, path = stack.pop()
            if len(path) == 5:
                if path[0] in neighbor_list[current_atom]:
                    for k in range(5):
                        p1, p2 = path[k], path[(k + 1) % 5]
                        pentagon_bonds.add(tuple(sorted((p1, p2))))
                continue
            for neighbor in neighbor_list[current_atom]:
                if len(path) > 1 and neighbor == path[-2]: continue
                if neighbor not in path:
                    stack.append((neighbor, path + [neighbor]))

    # Define colors and create the color property for bonds.
    color_6_5 = (1.0, 0.4, 0.4); color_6_6 = (0.4, 0.4, 1.0)
    bond_color_property = []
    for p1, p2 in bond_topology:
        bond_tuple = tuple(sorted((p1, p2)))
        bond_color_property.append(color_6_5 if bond_tuple in pentagon_bonds else color_6_6)
    
    # Create the 'Color' property on the bonds. OVITO will use this automatically.
    data.particles_.bonds_.create_property('Color', data=bond_color_property)
    
    print(f"Frame {frame}: Data coloring complete.")

# --- 1. Setup File Paths ---
# IMPORTANT: Update this path to point to your C60.xyz file.
xyz_file_path = '/Users/kannansekarannuradha/Downloads/C60.xyz'
output_dir = os.path.dirname(xyz_file_path)
image_path = os.path.join(output_dir, "c60_nanocluster_bonds.png")
movie_path = os.path.join(output_dir, "c60_rotation.mp4")

# --- 2. Create and Compute a Temporary Data Processing Pipeline ---
print("Creating temporary data processing pipeline...")
processing_pipeline = Pipeline(source=import_file(xyz_file_path))

# Add modifiers to the temporary pipeline
processing_pipeline.modifiers.append(CreateBondsModifier(cutoff=1.9))
processing_pipeline.modifiers.append(PythonScriptModifier(function=setup_coloring))

# Compute the temporary pipeline to get the final, colored data object.
print("Computing data pipeline...")
final_data = processing_pipeline.compute()
print("Pipeline computed.")

# --- 3. Create a New, Clean Pipeline for Rendering ---
# CORRECTED: To prevent recursion errors, we create a new, clean pipeline
# for rendering, using the data we just computed as a static source.
# This pipeline has no modifiers and is safe to render.
print("Creating clean pipeline for rendering...")
render_pipeline = Pipeline(source=StaticSource(data=final_data))
render_pipeline.add_to_scene() # Add the clean pipeline to the scene.

# --- 4. Setup Visual Appearance ---
# Set visual properties on the data object within the clean pipeline.
print("Setting visual properties...")
render_pipeline.source.data.particles.vis.radius = 0.35
render_pipeline.source.data.particles.bonds.vis.width = 0.2

# --- 5. Setup Viewport and Camera ---
vp = Viewport()
vp.type = Viewport.Type.Perspective
# The scene now contains our clean render_pipeline.
vp.scene = render_pipeline # Explicitly tell the viewport to look at our clean scene
vp.camera_pos = (12, 12, 12)
vp.camera_dir = (-1, -1, -1)
vp.fov = np.deg2rad(30.0)
vp.zoom_all()

# --- 6. Render Static Image ---
print(f"Rendering image to: {image_path}")
vp.render_image(
    filename=image_path,
    size=(1200, 1200),
    renderer=TachyonRenderer()
)
print("Image rendering complete.")

# --- 7. Render Rotation Movie ---
print(f"Rendering movie to: {movie_path}")

# Define a Python generator function to animate the camera for the movie renderer.
initial_pos = vp.camera_pos
total_frames = 120
def animate_camera(frame: int, total_frames: int):
    # Calculate rotation angle for the current frame
    angle = 2 * np.pi * frame / total_frames
    
    # Rotate camera position around the Z-axis
    x = initial_pos[0] * np.cos(angle) - initial_pos[1] * np.sin(angle)
    y = initial_pos[0] * np.sin(angle) + initial_pos[1] * np.cos(angle)
    z = initial_pos[2]
    vp.camera_pos = (x, y, z)
    
    # Report progress to the console
    if (frame+1) % 15 == 0:
        print(f"  ...rendered movie frame {frame+1}/{total_frames}")

    # Yield control back to the renderer for this frame.
    yield

# Use the render_movie() function with the animation generator.
vp.render_movie(
    movie_path,
    total_frames,
    animate_camera,
    fps=30,
    size=(1080, 1080),
    renderer=TachyonRenderer()
)

print("Movie rendering complete.")
print("\nScript finished successfully.")
