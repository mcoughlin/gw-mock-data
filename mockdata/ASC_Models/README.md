# ASC Models

Aim is to provide realistic ASC control signal models so that the designed neurel network architecture is accurate enough to carry out the real-data ready bilinear noise subtraction.


## Reading the model 
L1-ASC_DHARD_P.npz added on  [20/7/2017] 

* Actual Spectrum : https://ldvw.ligo.caltech.edu/ldvw/view?act=getImg&imgId=173166
* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/ASC_Modeled_Response.png


Can be read in python in the following manner:
```
asc_model = np.load(ASC_model.npz)
zeros = asc_model['z']
poles = asc_model['p']
gain  = asc_model['k']
```


H1-ASC_CHARD_P.npz added on 21/7/2017 

* Actual Spectrum : https://ldvw.ligo.caltech.edu/ldvw/view?act=getImg&imgId=173282
* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/H1_CHARD_Model.png  

# Lower order ASC Models 

These ones are mimic response from a smoothed version of the ASC control signal spectra. All major resonances are removed to reduce the model order. Unlike the previous  models which were manually built using MATLAB constolSystemDesigner tool, these ones are generated via an automatic scripts and hence only approximate the actual response.
All the models use a maximum of 12 poles to approximate the response. Most added models fit the data well above 1 Hz.

To generate your own model for the measured data, run ASC_Control_Signal_Modeling.m ( ZPK .mat files are converted to .npz using convert_zpk_matlab_python.py)
 
L1-ASC-DHARD_Y_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-DHARD_Y_PeakSmoothed.png

L1-ASC-DHARD_P_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-DHARD_P_PeakSmoothed.png

L1-ASC-CHARD_Y_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-CHARD_Y_PeakSmoothed.png

L1-ASC-CHARD_P_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-CHARD_P_PeakSmoothed.png
 

H1-ASC-DHARD_Y_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-DHARD_Y_PeakSmoothed.png
 

H1-ASC-DHARD_P_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-DHARD_P_PeakSmoothed.png

H1-ASC-CHARD_Y_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-CHARD_Y_PeakSmoothed.png 

H1-ASC-CHARD_P_PeakSmoothed.npz 

* Modeled Spectrum : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/L1-ASC-CHARD_P_PeakSmoothed.png 

