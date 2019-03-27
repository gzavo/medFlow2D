# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 09:26:11 2015

@author: Gabor Zavodszky
"""

from matplotlib.pyplot import imshow 
import numpy as np
from PIL import Image, ImageDraw




### Set these! ###
###  Use [mm]! ###
imgRes = 512         # Resolution of the final image
D = 2.0              # Diameter of the vessel
L = 1.5              # Arc-width of the aneurysmal neck
r = 1.5              # Radius of the aneurysmal sac
w = 2.0              # Width of the stent in pixels
Length = 8			 # Number of diameters for parent artery length
oFile = 'geom2.ppm'   # Name of the output file
pFile = 'param.ini'  # Parameters of this geometry
##################



# Calculate the maximum extents of the image

imgMinWidth = Length*D
dx = imgMinWidth / (imgRes-1.0) 
#imgMinHeight = (D + r * (1.0 + np.sin(np.arccos(0.5*L/r))))

print "Pixel resolution for setting up the simulation [m]: ", (dx / 1000.0)

# Set parameters in pixel units
pD = D / dx
pL = L / dx
pr = r / dx
pwidth = imgMinWidth / dx
#pheight = imgMinHeight / dx


# Creating the image
img = Image.new( 'RGB', (imgRes,imgRes), "black") # create a new black image
pixels = img.load() # create the pixel map


# Draw vessel
for i in range(1, int(pwidth)):    # for every pixel:
    for j in range(1, int(pD)):
        if(j<pD):
            pixels[i,j] = (255, 255, 255)

# Calc aneurysm neck coordinates
p1 = pD
p2 = (0.5*(pwidth-pL))

# Fit circle to two points and a radius
x1 = p1
y1 = p2
x2 = p1
y2 = p2 + pL

# Center of the circle
cx1 = (y1 + y2)/2.0
cy1 = x1 + np.sqrt((pr**2) - (0.5*pL)**2)

# Draw aneurysm
for i in range(img.size[0]):    # for every pixel:
    for j in range(img.size[1]):
        if( ((i-cx1)**2 + (j-cy1)**2) < (pr)**2 ):
            pixels[i,j] = (255, 255, 255)

# Draw inlet and outlet
for i in range(1, int(np.floor(pD))):
    pixels[pwidth ,i] = (255, 0, 0)

for i in range(1, int(np.floor(pD))):
    pixels[0, i] = (0, 255, 0)

# Draw stent
for i in range(int(y1-2*w), int(y2+1+2*w)):
    if (int(w)%2 == 0):
        for j in range(int(pD), int(pD+w)):
            if (pixels[i, j] == (255, 255, 255)):
                pixels[i, j] = (0, 0, 255)

    else:
        for j in range(int(pD), int(pD+w)):
            if (pixels[i, j] == (255, 255, 255)):
                pixels[i, j] = (0, 0, 255)


#draw = ImageDraw.Draw(img)
#draw.line([(y1, x1),(y2, x2)], fill=(0,0,255), width=6)

# Rotate, crop and save image
out = np.array(img.transpose(Image.ROTATE_180)).astype(np.uint8)

#imshow(out)
ifile = Image.fromarray(np.array(out[(imgRes-3-int(cy1+pr)):]).astype(np.uint8), mode='RGB')
ifile.save(oFile)

# Write out parameters in [m]
with open(pFile, 'w') as f:
    f.write('[geometry]\n')
    f.write('dx=%e\n' % (dx/1000))
    f.write('imgResolution=%d\n' % (imgRes))
    f.write('vesselDiameter=%e\n' % (D/1000))
    f.write('aneurysmNeck=%e\n' % (L/1000))
    f.write('aneurysmR=%e\n' % (r/1000))