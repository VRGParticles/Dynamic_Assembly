# Dynamic Assembly
This is a repository for the three custom algorithms used in: 

[*Tools for the quantification of the dynamic assembly of colloidal chains of ellipsoidal particles*](https://doi.org/10.1016/j.colcom.2022.100661)

This repository provides image analysis for videos of polymorphic particle/colloid assemblies.  
Particles are classified into groups based on chain counts. Particles are then passed to TrackPy to track the particles over time. The classification and tracking data is then combined to study particle dynamics.

# What is Included #

**MATLAB Files (.m):**

*Parameters_Particle_Radar.m* Parameters to edit to run the main script.

*Particle_Radar.m* The main script.

*find_angle_differences.m* Finds the angle differences between particles.

*find_parallel_neighbors.m* Finds parallel particle neighbors.

*image_preprocess.m* Preprocesses images from raw microscope image to binary.

*process_tracking_data.m* Process tracking data from TrackPy.

*Particle_Radar_Trackpy.ipynb* Jupyter notebook used to run TrackPy.


### General Algorithm Overview

## Getting Started
*MATLAB Requirements*

To use these files, you must have MATLAB installed.  Scripts were written using MATLAB R2021a. These algorithms use various toolboxes that must also be installed.

*Other Prerequisites*


[Trackpy](http://soft-matter.github.io/trackpy/v0.5.0/)

and 

[export_fig](https://github.com/altmany/export_fig/releases/tag/v3.21)


*Other Input Files*

The raw microscope video was converted to a stack of frames. 

## Running the Scripts

All input variables are located in the parameters file.



### More Advanced Script Customization 
For more advanced image analysis,  i.e. customizing these algorithms more toward a user-provided image instead of the provided input files, the user can edit the script. It is imagined that the type of particle assembly formation analyzed can be changed (i.e., particles form a different structure than chains).


### Troubleshooting
MATLAB related issues can be solved by viewing the [MATLAB documentation.](https://www.mathworks.com/help/index.html)

### Author

* **Veronica Grebe** - [Weck Group, New York University](http://weckresearchlab.com)

For a full list of manuscript authors who contributed to the overall project, see the [manuscript](https://doi.org/10.1016/j.colcom.2022.100661).

### License

This project is licensed under the GNU GPL-3.0 License - see the [LICENSE](LICENSE) file for details
