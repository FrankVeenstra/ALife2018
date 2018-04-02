# ALife2018
Minimal working example of a genetic algorithm on spatial and non-spatial model using H-IFF as the optimization function. 

This is a temporary version of the code used to run an evolutionary algorithm on a spatial model and a non spatial model. 
The code in its current form contains many messy parts but it will be updated to a more decent version for a camera ready paper. 

The repository contains two folders, the GAbit folder and the Cellular War folder. 
In the GAbit folder you can find a manual implementation of a genetic algorithm using the H-IFF function. 
Simply open GAbit/GAbit.py to do an example run. 
The uploaded code is a minimal version, other types of algorithms have also been implemented but are out of the scope of the paper. 
When running the python file it will commence 4 evolutionary runs for four different types of mortality rates: 0.0,0.02,0.04,0.06,0.08 and 0.1 on 32 bit H-IFF to show the performance difference. 
These runs take approximately 5 minutes to complete. 

The SpatialModel folder contains the spatial model which is adapted from an old version of Cellular war uploaded by Muzkaw: https://github.com/Muzkaw/Cellular-War/.
To run the spatial model program simply open "SpatialModel/Release/SpatialModel.exe". 
The settings.csv file can be altered to change the settings of the sptatial model. 
The data of the spatial model will be stored in another .csv file which is called test.csv by default. 
The spatial model uses the Simple and Fast Multimedia Library (SFML) which can be downloaded here: https://www.sfml-dev.org/ .

Contact me for questions or to get access to extended versions of the programs.  
