import numpy as np
import matplotlib.pyplot as plt

np.random.seed(999999999)
image = np.random.normal(1300, 10, (5,5))
# creating 5 x 5 array with mean and std
image[2,2] += 180
#selecting centre pixel and adding 180 to create star

image = image - np.median(image)
# removing bkdg of image


plt.imshow(image, cmap= "gray_r", interpolation= "none")
plt.show()
#shows image
star = image[2,2]
#selecting star
bkdg = np.delete(image, 12)
#deleting star

std = np.std(bkdg)
# finding std of bkdg aka finding noise

SNR = star / (4 * std)
#calculates signal to noise ratio
print(SNR)

images = np.zeros((5,5))
#defining images and making an array of zeros

for i in range (10):
#creating for loop
    image = np.random.normal(1300, 10, (5,5))

    image[2,2] += 180

    images += image
#to keep creating image (copy from above)

images = images - np.median(images)




star = images[2,2]
bkdg = np.delete(images, 12)


std = np.std(bkdg)


SNR = star / (4 * std)
print(SNR)

######
images = np.zeros((5,5))
#defining images and making an array of zeros
SNRlist = []

for i in range (500):
#creating for loop
    image = np.random.normal(1300, 10, (5,5))

    image[2,2] += 25.3

    images += image
#to keep creating images (copy from above)

    images = images - np.median(images)




    star = images[2,2]
    bkdg = np.delete(images, 12)


    std = np.std(bkdg)


    SNR = star / (4 * std)
    SNRlist.append(SNR)
xlist = np.arange(1,501, 1)
a = 0.7
plt.plot(xlist, SNRlist)
f = a * np.sqrt(xlist)
plt.plot(xlist, f, label = a)
plt.xlabel("number of exposures")
plt.ylabel("SNR")
plt.title("SNR for 21 mag star at LFOP")

plt.legend()
plt.savefig("statslab.png")
plt.show()
