import cv2
import numpy as np
import json
import random
import os
import argparse
import mitsuba
mitsuba.set_variant('scalar_rgb')
from mitsuba.core import Thread
from mitsuba.core.xml import load_file

ap = argparse.ArgumentParser()
ap.add_argument("-p", "--path", required=True,
                help="Path to .xml file")
ap.add_argument("-o", "--output", required=True,
                help="Path to output file")
ap.add_argument("-spp","--spp",required=True,
                help="Number of samples per pixel")

args = vars(ap.parse_args())

scene_path = args['path']
output_path = args['output']
spp_p = args['spp']

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(scene_path))

scene = load_file(scene_path, spp = spp_p)
# Get the scene's sensor (if many, can pick one by specifying the index)
sensor = scene.sensors()[0]

# Call the scene's integrator to render the loaded scene with the desired sensor
scene.integrator().render(scene, sensor)
# The rendered data is stored in the film
film = sensor.film()

# Write out a tone-mapped JPG of the same rendering
from mitsuba.core import Bitmap, Struct
img = film.bitmap(raw=True).convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
img.write(os.path.join(output_path,f'{spp_p}spp_mitsuba.png'))