# OVITO Pro Python script for C60 nanocluster visualization
# -----------------------------------------------------------
#
# *** IMPORTANT ***
# This script is formatted as a VIEWPORT OVERLAY.
#
# *** HOW TO USE ***
# 1. FIRST, manually load your 'C60.xyz' file into OVITO.
# 2. Go to the "Overlays" tab in the command panel (on the right).
# 3. Click "Python script" to add a new script overlay.
# 4. Load this Python script file.
# 5. The script will automatically run and color the bonds of your loaded file.

import ovito
from ovito.modifiers import CreateBondsModifier, PythonScriptModifier
import numpy as np

# --- Python Script Function for Direct Data Coloring ---
# This is the core logic that will be embedded in the PythonScriptModifier.
def color_bonds_by_ring_type(frame, data):
    """
    This function finds rings and directly creates 'Color' properties
    for both particles and bonds.
    """
    # This check is important because the modifier might run before bonds exist.
    if not data.particles.bonds:
        return

    # Set a uniform color for all particles.
    num_particles = data.particles.count
    particle_color = (0.8, 0.8, 0.8) # Neutral grey
    data.particles_.create_property('Color', data=[particle_color] * num_particles)

    # Access the bond topology directly.
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
    
    data.particles_.bonds_.create_property('Color', data=bond_color_property)
    
    # After coloring, we also set the visual appearance from within the modifier.
    data.particles.vis.radius = 0.35
    data.particles.bonds.vis.width = 0.2


# --- Main Overlay Function ---
# This function is called by OVITO when rendering the viewport.
# The context parameter must have a default value for compatibility.
def render(painter, context=None):
    if context is None: return

    try:
        pipeline = context.scene.pipelines[0]
    except IndexError:
        return # No pipeline loaded, do nothing.

    # Check if our specific Python script modifier has already been added.
    # This is the safest way to prevent the script from running multiple times.
    modifier_exists = any(
        isinstance(m, PythonScriptModifier) and m.function.__name__ == 'color_bonds_by_ring_type' 
        for m in pipeline.modifiers
    )

    if not modifier_exists:
        print("Overlay script: Adding modifiers for bond coloring...")
        
        # We simply add the modifiers to the pipeline.
        # OVITO's main application will handle computing them.
        # This avoids any forbidden operations within the overlay script.
        pipeline.modifiers.append(CreateBondsModifier(cutoff=1.9))
        pipeline.modifiers.append(PythonScriptModifier(function=color_bonds_by_ring_type))

        # We ask the viewport to redraw to make sure the changes are picked up.
        if context.viewports:
            context.viewports[0].update()

        print("Overlay script: Modifiers added. OVITO will now process the data.")

