#!/usr/bin/python3

from PIL import Image, ImageOps, ImageDraw, ImageFilter
import numpy as np
from matplotlib import pylab as plt

im1 = Image.open('bs_qe.png')

im2 = np.array(Image.open('bs_dftb.png'))
im2_swap = im2.copy()
im2_swap[:,:,1] = im2_swap[:,:,2] # blue -> green
im2_swap[:,:,2] = im2_swap[:,:,0] # red -> blue
im2_swap[:,:,0] = im2[:,:,1] # green -> red

Image.fromarray(im2_swap).save('bs_dftb_red.png')

im2v = Image.open('bs_dftb_red.png')
mask = Image.open('bs_qe.png').convert('1')

im1.paste(im2v,(0,0),mask)
im1.save('comparing.png')
