from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdu = fits.open("combined_lum.fit") #name of galaxy in parentheses
data = hdu[].data #index in brackets

plt.imshow(data,cmap = 'gray_r' , vmin = , vmax = )
plt.colorbar()
plt.show()

plt.hist(data.flat,50,(0,5000))
plt.show()

x = #coordinates of center of star you are looking for x
y = #coordinate y

x_background_min = #x coordinate at bottom right of graph
x_background_max = #xcoordinate at bottom left of graph +1
y_background_max = #y coordinate at top right of graph +1
y_background_min = #y coordinate at bottom right of graph

radius = #number of pixels from edge of star to center

star = data[y-radius:y+radius+1, x-radius:x+radius+1] #plus one because we otherwise would not get the last pixel
background = data[y_background_min:y_background_max,x_background_min:x_background_max]

bg = np.mean(background)

plt.imshow(star)
plt.colorbar()
plt.show()


flux = np.sum(star)

magnitude = -2.5*np.log10(flux) + offset

print(flux,magnitude)
