# Compare scanning strategies

This Python script takes as input several sets of Time Ordered
Information (TODs), as the ones produced by scientific spacecraft and
sky survey experiments, and it produces an animation which displays
the coverage of the sky as a function of time.

## Installation

The program requires Python 3 and the following libraries:

- click, to parse the command line
- matplotlib, for low-level graphics stuff
- numpy, for number crunching
- healpy, to project the pointings into Healpix maps
- astropy, to load the FITS files containing the TODs

## Usage

The program requires an INI file as its only input. The output of the
program is a set of images that can be built together to create a
video.

## How to make a MP4 file?

If you have `ffmpeg`, you can use the following command to combine the
files `img0000.png`, `img0001.png`, etc. into the file `animation.mp4`:

    ffmpeg -i img%04d.png -c:v libx264 -pix_fmt yuv420p -r 24 animation.mp4

In this example, the number of frames per second (FPS) is set to 24,
which is a fairly common standard.

## License

The program is distributed under a MIT license. See the file
`LICENSE`.
