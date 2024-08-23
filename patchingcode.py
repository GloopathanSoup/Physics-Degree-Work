# -*- coding: utf-8 -*-
"""PatchingCode.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1XoSBtPDPDV8eet1uysg4hooDdcwbdfk0
"""

#Assigning values for testing

BigImageRGBlocation = '1.jpeg'

BigImageMaskLocation = 'mask1.png'

#InputtedBigImageLocation = 'Unico.png'

#ExportLocation = '.'

#os.getcwd()

#1 is iceberg 0.8 is icesheet 0.6 is land

!zip -r /content/patchexport2.zip /content/DATA2

from google.colab import files
files.download("/content/patchexport2.zip")

!pip install Patchify
!pip install pillow

import patchify as pat
import os
import numpy as np
import matplotlib.pylot as plt
import random
from PIL import Image

#TRAINING DATA PREPROCESSING

#This is intended to be run once per satellite scene. For industrial implementation, this should be modified into a function for scalability, though it suffices for small numbers of scenes as was the case in this project.

#Loads the big images
BigImageRGB = Image.open(BigImageRGBlocation).convert('RGB')
BigImageMask = Image.open(BigImageMaskLocation).convert('RGB')

#splits the big images into many smaller 500x500 images ('patches')
RGBpatches = pat.patchify(np.asarray(BigImageRGB), (500, 500, 3), step = 500)
MaskPatches = pat.patchify(np.asarray(BigImageMask), (500, 500, 3), step = 500)

#Exporting patches

c = 0 #counter to split into train test and validation

#RGB patch exporting
for i in range(RGBpatches.shape[0]): #Loops over every patch of the inputted image
  for j in range(RGBpatches.shape[1]):
    if c%10 == 0:
      filename = "/content/DATA2/test/RGBtest/RGBimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(RGBpatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1
    elif (c + 1)%10 == 0:
      filename = "/content/DATA2/val/RGBval/RGBimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(RGBpatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1
    else:
      filename = "/content/DATA2/train/RGBtrain/RGBimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(RGBpatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1

c = 0

#Mask patch exporting
for i in range(MaskPatches.shape[0]): #Loops over every patch of the inputted image
  for j in range(MaskPatches.shape[1]):
    if c%10 == 0:
      filename = "/content/DATA2/test/masktest/maskimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(MaskPatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1
    elif (c + 1)%10 == 0:
      filename = "/content/DATA2/val/maskval/maskimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(MaskPatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1
    else:
      filename = "/content/DATA2/train/masktrain/maskimage1_{}_{}.jpg" #IMPORTANT: change image1 to image2 etc for every large image processed
      im = Image.fromarray(MaskPatches[i, j, 0]).convert('RGB')
      im.save(filename.format(i, j))
      c += 1

#PREDICTION PROCESSING

#loads the new image we are segmenting
InputtedBigImage = Image.open(InputtedBigImageLocation).convert('RGB')

#slices the new image into patches
InputtedImagePatches = pat.patchify(np.asarray(InputtedBigImage), (500, 500, 3), step = 500)

#creating a tuple recording the size and number of channels of the original image for use in unpatching the masks
sizetuple = InputtedBigImage.size
sizelist = list(InputtedBigImage.size)
sizelist.append(3) #adding the channel dimension
TrueSizeTuple = tuple(sizelist)

#Generating and storing new masks for every patch
GeneratedMaskPatches = np.empty_like(InputtedImagePatches) #Creates an object to store the generated masks of every patch
for i in range(InputtedImagePatches.shape[0]): #Loops over every patch of the inputted image
  for j in range(InputtedImagePatches.shape[1]):
    GeneratedMaskPatches[i, j] = PredictionCode(InputtedImagePatches[i, j, 0, :, :, :]) #Uses your prediction function. Either that or put the entirety of that code in here

#Putting together our generated masks into one big mask for the whole scene
FinalGeneratedMask = pat.unpatchify(GeneratedMaskPatches, TrueSizeTuple)

#Visualisation
fig = plt.figure(figsize=(10, 7))

fig.add_subplot(1, 2, 1)
plt.imshow(InputtedBigImage)
plt.axis('off')
plt.title("Original Image")

fig.add_subplot(1, 2, 2)
plt.imshow(FinalGeneratedMask)
plt.axis('off')
plt.title("Segmentation Mask")

#Exporting
filename = "/content/finalmaskexport/generatedmask_{}.jpg"
im = Image.fromarray(RGBpatches[i, j, 0]).convert('RGB')
im.save(filename.format(random.choice(range(9999))))

#testing stuff, ignore

#loads the new image we are segmenting
InputtedBigImage = Image.open(InputtedBigImageLocation).convert('RGB')

#slices the new image into patches
InputtedImagePatches = pat.patchify(np.asarray(InputtedBigImage), (500, 500, 3), step = 500)

#creating a tuple recording the size and number of channels of the original image for use in unpatching the masks
sizetuple = InputtedBigImage.size
sizelist = list(InputtedBigImage.size)
sizelist.append(3) #adding the channel dimension
TrueSizeTuple = tuple(sizelist)

GeneratedMaskPatches = np.empty_like(InputtedImagePatches) #Creates an object to store the generated masks of every patch
for i in range(InputtedImagePatches.shape[0]): #Loops over every patch of the inputted image
  for j in range(InputtedImagePatches.shape[1]):
    GeneratedMaskPatches[i, j] = InputtedImagePatches[i, j, 0, :, :, :] #Hamza is this correct? Uses our (currently non-existent) prediction code to generate a list full of generated masks for every patch of the inputted image

FinalGeneratedMask = pat.unpatchify(GeneratedMaskPatches, TrueSizeTuple)

fig = plt.figure(figsize=(10, 7))

fig.add_subplot(1, 2, 1)
plt.imshow(InputtedBigImage)
plt.axis('off')
plt.title("Original Image")

fig.add_subplot(1, 2, 2)
plt.imshow(FinalGeneratedMask)
plt.axis('off')
plt.title("Segmentation Mask")

print(InputtedBigImage.size)
print(GeneratedMaskPatches.shape)
#plt.imshow(InputtedBigImage)
#plt.imshow(GeneratedMaskPatches[1, 1, 0])

#Testing stuff

InputtedBigImage = Image.open(InputtedBigImageLocation).convert('RGB')
InputtedImagePatches = pat.patchify(np.asarray(InputtedBigImage), (250, 250, 3), step = 250)
print(InputtedImagePatches.shape)
os.getcwd()

#PatchListToPatchify(PatchifyToPatchList(InputtedImagePatches)).shape
#InputtedImagePatches.shape

#Segments = PatchifyToPatchList(InputtedImagePatches)

PatchList = []
for i in range(InputtedImagePatches.shape[0]): #loops over every patch horizontally and vertically
  for j in range(InputtedImagePatches.shape[1]):
    GetPatch = InputtedImagePatches[i, j, 0, :, :, :] #collects the patch as a single image
    PatchList.append(GetPatch)

len(PatchList)

#plt.imshow(Segments[5])