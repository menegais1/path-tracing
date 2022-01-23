import cv2
import numpy as np
import json
import random
import os
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-im1", "--image1", required=True,
                help="Image 1 path")
ap.add_argument("-im2", "--image2", required=True,
                help="Image 2 path")
ap.add_argument("-fm", "--finalimage", required=True,
                help="Final Image path")
args = vars(ap.parse_args())

image1_path = args['image1']
image2_path = args['image2']
final_image_path = args['finalimage']

image1 = cv2.imread(image1_path)
image2 = cv2.imread(image2_path)

image1_mean = np.mean(image1.astype(np.float32), axis=2)
image2_mean = np.mean(image2.astype(np.float32), axis=2)
print(image1_mean.shape)
final_image = np.abs(image2_mean - image1_mean).astype(np.uint8)
final_image = np.stack((final_image,)*3, axis=-1)

final_image[:,:,0] = 0;
final_image[:,:,2] = 0;

cv2.imwrite(final_image_path, final_image.astype(np.uint8))

