# Mock Data

This project is for generating mock data to be used for testing our nonlinear regression tools.

**Please put your noise subtraction code into the NonlinearRegression area, not in this project.**

## The Concept
We think that there is a background noise level set by thermal and shot noise.
The dominant noise which limits LIGO in the 40-200 Hz band at the moment is not due to
these more fundamental noise sources, but is some other mysterious noise.

In order to see if this noise is due to a nonlinear combination of some measured
observeables in the system, we want to produce some mock data using combinations
of PEM, SUS, SEI, LSC, and ASC channels. We add this onto the Gaussian noise background
due to the thermal and shot noise.

Then we use a supervised learning algorithm to predict that noise based on the
input witness channels. If our algorithm is good it will be able to figure out how
to subtract the nonlinear noise from the interferometer.

## Code Organization

The mockdata package contains tools for generating mock gaussian background noise
(contained in mockdata/mock_bg.py) and nonlinear noise based on random witnesses
(mockdata/mock_noise.py). To generate and save example timeseries data, run makeNoise.py
or call mockdata.starting_data() within a script. For background only, call
mockdata.mock_bg.bucket_noise(). To generate complex spectrogram data,
run make_STFT_data.py. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Example

kinit albert.einstein@LIGO.org
python getDARMbilinearData.py L1

### Prerequisities

1. Python
1. Scipy
1. Scikit-learn
1. TensorFlow
1. pynoisesub (from https://git.ligo.org/NoiseCancellation/pynoisesub)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Claude Shannon
* Norbert Weiner
* LIGO
* LIGO Science Collaboration
* LIGO Virgo Collaboration
* Widrow
