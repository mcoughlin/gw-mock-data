# Quad SUS Models

* Quad Suspension ZPK Models generated using generate_QUAD_Model_Production.m. Currently has the following models : 
* Penultimate Stage Pitch drive [Nm] to Test Mass Pitch disp [in rad]
* Penultimate Stage Yaw drive [Nm] to Test Mass Yaw disp [in rad]

## Bode Plot
* Pitch to Pitch : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/SUS_model_PUM_Pitch_to_TST_Pitch.png 

* Yaw to Yaw     : https://ldas-jobs.ligo.caltech.edu/~nikhil.mukund/BiLinear/SUS_model_PUM_Yaw_to_TST_Yaw.png 


## Reading the model 
Can be read in python in the following manner:
```
model = np.load('SUS_model_PUM_Pitch_to_TST_Pitch.npz')
zeros = model['z']
poles = model['p']
gain  = model['k']
```





