#!/usr/bin/env python3
import math
import subprocess
from ovito.io import import_file
from ovito.modifiers import RingFinderModifier
from ovito.vis import (
    Viewport, TachyonRenderer, ColorCodingRenderingSettings
)

# 1) Load the C60 file
pipeline = import_file('/mnt/data/C60.xyz')

# 2) (Implicit) — your cluster is now in `pipeline`

# 3) Color by ring type (pentagons vs hexagons)
ring_mod = RingFinderModifier(ring_sizes=[5,6])
pipeline.modifiers.append(ring_mod)

color_settings = ColorCodingRenderingSettings()
color_settings.property       = 'Ring sizes'
color_settings.coloring_mode  = ColorCodingRenderingSettings.ColoringMode.Group
color_settings.use_colormap   = False
color_settings.group_color_map = {
    5: (1.0, 0.0, 0.0),  # pentagons → red
    6: (0.0, 1.0, 0.0),  # hexagons → green
}

# 4) Render a nice still image
vp = Viewport(type=Viewport.Type.Perspective)
vp.source      = pipeline
vp.renderer    = TachyonRenderer()            # high-quality Tachyon engine
vp.camera_pos  = (0.0, 0.0, 20.0)              # pull back a bit
vp.camera_dir  = (0.0, 0.0, -1.0)
vp.background_color = (1.0, 1.0, 1.0)          # white background

vp.render(
    filename       = 'C60.png',
    size           = (1200, 900),
    background     = vp.background_color,
    color_settings = color_settings
)

# 5) Produce a rotating movie
num_frames = 120
radius     = 20.0

for i in range(num_frames):
    angle = 2*math.pi * i / num_frames
    x = radius * math.cos(angle)
    z = radius * math.sin(angle)
    vp.camera_pos = (x, 0.0, z)
    # always look toward cluster center at origin:
    vp.camera_dir = (-x, 0.0, -z)
    vp.render(
        filename       = f'frame_{i:03d}.png',
        size           = (800, 600),
        background     = vp.background_color,
        color_settings = color_settings
    )

# Stitch frames into an MP4 (requires ffmpeg on your PATH)
subprocess.call([
    'ffmpeg', '-y',
    '-framerate', '30',
    '-i', 'frame_%03d.png',
    '-c:v', 'libx264',
    '-pix_fmt', 'yuv420p',
    'C60_rotation.mp4'
])

print("Done: still image → C60.png; movie → C60_rotation.mp4")
