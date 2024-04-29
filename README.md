# EphysUtils
**Work in progress.**

Pipeline to preprocess extracellular electrophysiology Neuropixels data acquired using SpikeGLX. 

### Notes about this pipeline:
- Works with the suite of SpikeGLX tools i.e. **CatGT**, **TPrime**, etc. : https://billkarsh.github.io/SpikeGLX/
- Borrows/adapts some code found in the SpikeGLX-adapted fork of the Allen's ecephys pipeline: https://github.com/jenniferColonell/ecephys_spike_sorting (e.g. mean waveform calculation)
- Written for Kilosort spike sorting (KS2 mostly):  https://github.com/MouseLand/Kilosort?tab=readme-ov-file (KS4 now)
- Includes semi-automated curation using Bombcell: https://github.com/Julie-Fabre/bombcell

### Overview of the pipeline
```mermaid
graph LR
    Step1["Event extraction <br/> filtering"] -.->|Next| Step2["Coil artifact <br/> correction"]
    Step2 -.->|Next| Step3(["Optional: <br/> chunk zeroing"])
    Step3 -.->|Next| Step4["Spike sorting <br/> & quality metrics"]
    Step4 -.->|Next| Step5["Data stream <br/> synchronization"]
    Step5 -.->|Next| Step6["Mean waveform <br/> & metrics"]
    Step6 -.->|Next| Step7["LFP <br/> analysis"]

````
#### Summary of the main steps

- **Events extraction (CatGT)**: extracts times of TTL pulses acquired with the NI card in the `nidq.bin` output file of SpikeGLX
- **Filtering (CatGT)**: common median referencing by default
- **Coil artifact correction (TPrime)**:
  1. synchronize extracted coil/whisker stimulation times to each IMEC probe base time
  2. at each artifact time, replace duration of artifact (3ms default) by mean voltage just before, for all channels
  3. create copy of .ap/.meta file with the "corrected" suffix 
- **Chunk zeroing (OverStrike)**: zero-out entire chunks of data in the recordings when there is unsalvageable noise
- **Spike sorting (Kilosort)**: spike sorting algorithm for neuron identification, calls Kilosort 2.0 from the Python MATLAB engine (see below)
    - Notes about spike sorting
- **Quality metrics**: runs quality metrics pipeline from **Bombcell** (CortexLab) from the MATLAB engine, with by modified default:
  - Plotting is set to off (one plot/cluster generated), set to True for initial debugging/inspection
  - Further splitting of non-somatic to mua/good is set False
  - Computations of ephys properties is set to False (not immediately necessary)
- **Data stream synchronization**
- 
#### Installation
- Install the associated conda environment:
- MATLAB e.g. R2021b - specify the MATLAB version to use when calling the MATLAB engine in Python:
  - In MATLAB command window, type `matlabroot` to get root path
  - In terminal, go to `<matlabroot>\extern\engines\pyton`, then type `python setup.py install`
  - If the previous did not work, try: https://ch.mathworks.com/matlabcentral/answers/1998578-invalid-version-r2021-when-installing-for-python-3-7-3-9
    
### Future improvements (and ideas):
- Adaptation/robustness for Neuropixels 2.0 probes specifications and metadata (all SpikeGLX & cie tools updated for each hardware changes) 
- Kilosort 4.0 called from python directly
- Integration of [SpikeInterface](https://github.com/SpikeInterface) tool(s)
- More LFP analyses...?

  
