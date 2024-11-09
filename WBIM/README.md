# WBIM
WBIM is an all-optical serial multiphoton microscope designed for sub-micrometer resolution 3D imaging of organ-scale complex tissue. This system combined multiphoton microscopy with femtosecond laser machining to iteratively image and section biological tissues of distinct optical and mechanical properties. We tested this system on in situ imaging of mouse brain-skull system. 

This repository contains a set of software and computational pipelines developed that enable multiday-long automated acquisition. Specifically: 
1. `Analysis` contains the the MATLAB script and classes for analyzing various experimental results. 
1. `Calibration` contains the MATLAB scripts and classes for instrument calibration, including: 
    1. Automated HWP-PBS-based laser power control calibration
    1. Automated Pockels-cell-PBS-based laser power control calibration 
    1. Servo motorized stage controler delay correction 
    1. Femtosecond laser machining characterization 
1. `Control` contains the major classes implemented for microscope control, including: 
    1. `WBIMConfig`: base class recording instrument parameters, coordinate transform, and configuration space. 
    1. `WBIMControlBase`: child class of `WBIMConfig`. Implements logging and motorized stage and piezo control. 
    1. `WBIMImaging`: child class of `WBIMControlBase`. Implements methods related to automated serial imaging, including: 
        1. Acquisition grid generation
        1. Imaging related hardware (laser, imaging HWP, imaging shutter) initialization and control
        1. Asynchronous data processing TCP server setup and data management
        1. ROI detection, classification, and acquisition trajectory path planning 
        1. Visual-based quality control and hardware abnormal state detection (scanner desynchronization, laser failures)
    1. `WBIMAblation`: child class of `WBIMControlBase`. Implements methods related to automated femtosecond laser maching, inculding: 
        1. Ablation related hardware (laser gating, ablation beam shutter, ablation HWP, pump, refractometers) initialization and control
        1. Ablation ROI detection (Multichannel 3D image spectrum unmixing, segmentation, classification), coverage path planning (TSP), ablation control sequence synthesis, execution and monitoring. 
    1. `WBIMControl`: child class of `WBIMControlBase`, integrates `WBIMImaging` and `WBIMAblation` for iterative imaging and ablation. 
    1. Additional classes include: 
        1. Instrument control: `ZaberController`, `WBIMPowerHWP`, `ThorlabsPowerMeter`
        1. Data management: `WBIMTileManager`, `WBIMImagingGrid2D`
        1. Email and message notification: `WBIMNotification`
        1. Ablation control related: `WBIMAblationPath3D`, `WBIMAblationPath2D`
1. `Laser_control` contains the classes for controlling and monitoring Coherent Chameleon laser 
1. `Postprocessing` contains the scripts and classes for image processing, asynchronous data processing, and data management. 
1. `Refractometer` contains the classes for controlling inline refractometer and adjusting solution index



