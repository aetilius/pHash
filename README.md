# pHash&trade; Perceptual Hashing Library

pHash is a collection of perceptual hashing algorithms for image,
audo, video and text media.  

## Installation

Build and install like so:

```
cd <project-dir>
cmake .
make
sudo make install
```

To run tests:

```
ctest
```
Look at the tests in the tests/ folder to see how they are used.


The build can exclude various hashing algorithms by passing
the variables USE_IMAGE_HASH, USE_AUDIO_HASH, USE_VIDEO_HASH
OR USE_TEXT_HASH to cmake program: `cmake -DUSE_IMAGE_HASH=1`

The Windows build currently does not work.  If anyone manages to
get it to work, please send a pull request with the changes.


## Image Hashing

There are four image hash algorithms in the pHash librarie.

1. DCT image hash\
   A 64-bit binary hash based on the global discrete cosine transform (dct)
   of the image.  Comparisons are made by calculating the hamming distance,
   or the number of places where the bits differ, between two hashes.

2. Radial Image Hash\
   A perceptual hash based on the variances of pixels taken from strips of
   various angles across the center of the image.
   
3. MH Image hash\
   Based on the mexican hast wavelet function.
   
4. BMB Image Hash - Box-Mean-Binarized Hash\
   The image is divided up into 16x16 pixel-sized boxes.  The
   mean value is calculated for each box in the gird, and a binarized
   hash is formed based on each mean value with respect to the median
   of those mean values. A basic hamming distance is used to compute
   distances.

There is an example program in `examples/imghash.cpp` to evaluate the
image hashes. Basically, it takes two directories of image files.  One
directory must be an exact replica of the other except for a specific
distortion applied to the images in one of the directories.  The
original/distorted pairs must have the same filename.

The program calculates the image hashes for each pair of images and
computes the distance.  A histogram of these distances is then plotted.

As a point of reference, a histogram of distances between the image
hashes of random dissimilar images is also plotted.  The peaks of the two histograms
should exhibit good separation. The program, imghash prints out the two histograms
to an output file.  This file can then be used by a tool such as gnuplot to plot
the histograms.

Here is a python script for processing a test set of image files:
[test set of image files](https://github.com/starkdg/pyConvnetPhash/blob/master/preprocess_image_files.py)


Run `./imghash --help` to show a list of options to provide.


### Results

Test results for the image hashes can be summed up in the following table. To conduct the test,
an assortment of about one hundred natural images was amassed.  The distorted replicas were then
created by the above script.  The entries in the table represent the percentage of similar distances
that fall above the threshold values.  The threshold values are set to be the average dissimilar distance
minus one and two standard deviations.  Values close to zero indicate the hashes ability to distinguish
between similar and wholely different images for that particular distortion.  

![image hash test result](https://github.com/starkdg/phash/raw/master/resources/imghash-test-results.png "image phash tests")

As you can see, the distortions proving the toughest are rotation, shear, crop, horizontal and vertical flip.  The dct image
hash offers fairly good results - especially given that it is only 64 bits of information.  The 32-byte length BMB hash offers
some modest improvements.  

The radial hash disappoints with a further degradation in the bright, dark and histogram equalization distortions.

## Audio Hash

An audio hash is included here for completeness.  There is a fuller java
implementation here: [JPHashAudio](https://github.com/starkdg/JPhashAudio)

## Video Hash

A video hash is included, although admittedly a bit buggy.  A fuller c++
implementation is available here: [ClipSeekr](https://github.com/starkdg/clipseekr)

## Text Hash




## Requirements for Example Programs

Boost 1.70.0 - Only the filesystem, system and program_options libraries.


## Requirements for Image Hash

CImg v.6.6 (included) \   
libtiff-dev (optional) \
libpng-dev (optional) \
libjpeg-dev (optional)


## Requirements for Video Hash

ffmpeg v2.8.15 "feynman release" \
  libavformat \
  libavcodec \
  libavswscale

## Requirements for Audio Hash

libsndfile v1.0.27 \
libsamplerate v0.1.8 \
libmpg123 v1.23.8 (optional) \
libvorbis v1.3.5 (optional) \
libogg v1.3.2 (optional)

