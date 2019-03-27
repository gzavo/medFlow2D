# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 09:26:11 2015

@author: Gabor Zavodszky
"""

from matplotlib.pyplot import imshow 
import numpy as np

# recent Pillow
from PIL import Image
import ImageDraw
# Alternatively
#from PIL import Image, ImageDraw


### Set these! ###
###  Use [mm]! ###
imgRes = 512         # Resolution of the final image
R = 3.5              # Mean radius of the vessel
D = 1.0              # Diameter of the vessel
L = 1.0              # Arc-width of the aneurysmal neck
r = 1.5              # Radius of the aneurysmal sac
oFile = 'geom.ppm'   # Name of the output file
pFile = 'param.ini'  # Parameters of this geometry
##################

def arc(draw, bbox, start, end, fill, width=1, segments=100):
    """
    Hack that looks similar to PIL's draw.arc(), but can specify a line width.
    """
    # radians
    #start *= np.pi / 180
    #end *= np.pi / 180
    
    # angle step
    da = (end - start) / segments

    # shift end points with half a segment angle
    start -= da / 2
    end -= da / 2

    # ellips radii
    rx = (bbox[2] - bbox[0]) / 2
    ry = (bbox[3] - bbox[1]) / 2

    # box centre
    cx = bbox[0] + rx
    cy = bbox[1] + ry

    # segment length
    l = (rx+ry) * da / 2.0

    for i in range(segments):

        # angle centre
        a = start + (i+0.5) * da

        # x,y centre
        x = cx + np.cos(a) * rx
        y = cy + np.sin(a) * ry

        # derivatives
        dx = -np.sin(a) * rx / (rx+ry)
        dy = np.cos(a) * ry / (rx+ry)

        draw.line([(x-dx*l,y-dy*l), (x+dx*l, y+dy*l)], fill=fill, width=width)

# Calculate the maximum extents of the image
toAneurysmCenter = (R + (D/2.0) + r*np.cos(np.arcsin(0.5*L/r))) * np.sin(np.pi/4)
imgMinWidth = toAneurysmCenter + r
dx = imgMinWidth / (imgRes-1.0) 

print "Pixel resolution for setting up the simulation [m]: ", (dx / 1000.0)

# Set parameters in pixel units
pR = R / dx
pD = D / dx
pL = L / dx
pr = r / dx

pRM = pR+pD/2.0
pRm = pR-pD/2.0

# Creating the image
img = Image.new( 'RGB', (imgRes,imgRes), "black") # create a new black image
pixels = img.load() # create the pixel map


# Draw vessel
for i in range(img.size[0]):    # for every pixel:
    for j in range(img.size[1]):
        if( (i**2 + j**2) < (pRM)**2 ):
            pixels[i,j] = (255, 255, 255)
        if( (i**2 + j**2) < (pRm)**2 ):
            pixels[i,j] = (0, 0, 0)


# Calc aneurysm neck coordinates
alpha = 0.5*pL/(pRM)
p1 = pRM * np.sin(np.pi/4.0 - alpha)
p2 = pRM * np.cos(np.pi/4.0 - alpha)

# Fit circle to two points and a radius
x1 = p1
y1 = p2
x2 = p2
y2 = p1

d = np.sqrt((x2-x1)**2 + (y2-y1)**2)
x3 = (x1+x2)/2.0
y3 = (y1+y2)/2.0
cx1 = x3 + np.sqrt(pr**2-(d/2)**2)*(y1-y2)/d
cy1 = y3 + np.sqrt(pr**2-(d/2)**2)*(x2-x1)/d	
#cx2 = x3 - np.sqrt(r**2-(d/2)**2)*(y1-y2)/d
#cy2 = y3 - np.sqrt(r**2-(d/2)**2)*(x2-x1)/d

# Draw aneurysm
for i in range(img.size[0]):    # for every pixel:
    for j in range(img.size[1]):
        if( ((i-cx1)**2 + (j-cy1)**2) < (pr)**2 ):
            pixels[i,j] = (255, 255, 255)


# Draw inlet and outlet
for i in range(int(np.floor(pRm+1.0)), int(np.floor(pRM)+1.0)):
    pixels[0,i] = (255, 0, 0)

for i in range(int(np.floor(pRm+1.0)), int(np.floor(pRM)+1.0)):
    pixels[i,0] = (0, 255, 0)


# Draw stent
draw = ImageDraw.Draw(img)
arc(draw, (-pRM-1, -pRM-1, pRM+1, pRM+1), np.pi/4 - (0.5*pL/pRM)-(1./pRM), np.pi/4 + (0.5*pL/pRM)+(2./pRM), (0,0,255), width=2, segments=101)

# Rotate and save image
out = np.array(img.transpose(Image.ROTATE_90)).astype(np.uint8)
#imshow(out)
ifile = Image.fromarray(np.array(out).astype(np.uint8), mode='RGB')
ifile.save(oFile)

# Write out parameters in [m]
with open(pFile, 'w') as f:
    f.write('[geometry]\n')
    f.write('dx=%e\n' % (dx/1000))
    f.write('imgResolution=%d\n' % (imgRes))
    f.write('vesselRadius=%e\n' % (R/1000))
    f.write('vesselDiameter=%e\n' % (D/1000))
    f.write('aneurysmNeck=%e\n' % (L/1000))
    f.write('aneurysmR=%e\n' % (r/1000))

