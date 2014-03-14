import numpy as np
import pyfits
import matplotlib as ml
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from PIL import Image, ImageDraw

""" This program analyzes the light distribution of M86 galaxy from its image by finding and plotting centroid and
    the principal axes """

def solution(file = 'M86.fits'):
    f = pyfits.open(file)   #Read the .fits file using the pyfits library which turns it into a 2DNumpy array
    image = f[0].data
    height = len(image)
    width = len(image[0])
    Mean = np.mean(image)
    Std = np.std(image)     
    Threshold = Mean + Std  # Threshold value to remove noise 
    for i in range(0, height):
        for j in range(0, width):
            if(image[i][j] < Threshold):
                image[i][j]=0
    x_centroid, y_centroid = findCentroid(image)    #Finding the centroid of the Light Distribution
    angle = findOrientation(image, x_centroid, y_centroid)  #Finding the angle of orientation
    plotAxes(image, angle, x_centroid, y_centroid)  # Plotting the results

""" Finds the centroid by performing weighted mean of the intensities for every pixel"""

def findCentroid(map):
    map_height = len(map)
    map_width = len(map[0])
    total_sum = np.sum(map)
    x_sum = 0
    y_sum = 0
    for i in range(0, map_height):
        for j in range(0, map_width):
            x_sum += j*map[i][j]
            y_sum += i*map[i][j]
    x_sum = x_sum / total_sum
    y_sum = y_sum / total_sum
    return x_sum, y_sum

""" Finds the angle of orientation by calculating second order moments """

def findOrientation(map, xcentroid, ycentroid):
    map_height = len(map)
    map_width = len(map[0])
    x_bar = xcentroid
    y_bar = ycentroid
    print "Centroid of the Light Distribution : " + str([x_bar, y_bar])
    mu_20 = 0
    mu_02 = 0
    mu_11 = 0
    for j in range(0, map_width):
        for i in range(0, map_height):
            mu_20+= ((j - x_bar)**2)*map[i][j]
    for j in range(0, map_width):
        for i in range(0, map_height):            
            mu_02+= ((i - y_bar)**2)*map[i][j]
    for j in range(0, map_width):
        for i in range(0, map_height):            
            mu_11+= ((j - x_bar)*(i - y_bar))*map[i][j]
    theta = 0.5*(np.arctan((2*mu_11)/(mu_20 - mu_02)))
    print "Orientation of the Light Distribution : " + str(theta)+ " radians"
    return theta                        

""" Plots the results with a log normalization to visualize important features  """    

def plotAxes(map, angle, x_centroid, y_centroid):
    hor = (x_centroid + 7*(math.cos(angle)))
    ver = (y_centroid - 7*(math.sin(angle)))
    hor1 =(x_centroid + 7*(math.cos(angle+ (math.pi)/2)))
    ver1 = (y_centroid - 7*(math.sin(angle+(math.pi)/2)))
    hor2 = (x_centroid + 7*(math.cos(angle+math.pi)))
    ver2 = (y_centroid - 7*(math.sin(angle+math.pi)))
    hor3 = (x_centroid + 7*(math.cos(angle+ math.pi+(math.pi)/2)))
    ver3 = (y_centroid - 7*(math.sin(angle+math.pi+(math.pi)/2)))
    imgplot = plt.imshow(map, cmap = cm.gray, origin = "lower"""",norm=lm(map.mean() + 0.5 * map.std(), map.max(), clip='True')""")
    plt.plot([hor2, hor],[ver2,ver], linestyle='-', linewidth=1, ms=5, mfc='blue', alpha=0.4)
    plt.plot([hor3, hor1],[ver3,ver1], linestyle='-', linewidth=1, ms=5, mfc='yellow', alpha=0.4)
    plt.plot([x_centroid],[y_centroid], 'o', ms=5, mfc='red', alpha=0.4)
    plt.title('GSoC 2014 solution : Krishna Chaitanya Chavati', color='black')
    plt.xlabel('X',color ='black')
    plt.ylabel('Y',color = 'black')
    plt.show()
    
    
if __name__ == '__main__':
    solution()
