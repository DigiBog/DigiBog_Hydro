# DigiBog_Hydro Model Documentation

## Overview

DigiBog_Hydro is a hydrological model that simulates water-table dynamics in peatland landscapes using the Boussinesq equation. The model is designed specifically for water movement simulation and does not include peat accumulation or loss processes.

### Key Features
- **Water-table dynamics simulation** over peatland landscapes
- **Management scenario analysis** (e.g. ditch drainage, ditch blocking)
- **2.5D finite difference approach** using the Boussinesq equation
- **Flexible boundary conditions** (Dirichlet, Neumann, inactive)
- **Variable peat properties** (hydraulic conductivity, porosity, thickness)
- **Climate forcing** via rainfall and potential evapotranspiration

### Model Specifications
- **Language:** Fortran 95
- **Units:** Centimetres and seconds
- **Licence:** GNU General Public License v3
- **Standards:** Modern Fortran

### Authors
- Paul J. Morris (p.j.morris@reading.ac.uk)
- Andy J. Baird (a.j.baird@leeds.ac.uk) - School of Geography, University of Leeds
- Lisa R. Belyea (l.belyea@qmul.ac.uk)
- Dylan M. Young (d.m.young@leeds.ac.uk)

For more information, visit: https://water.leeds.ac.uk/our-missions/mission-1/digibog/about-digibog/

## Basic Workflow

### 1. Model Setup
1. Define field geometry and grid dimensions
2. Set up boundary conditions and active cells
3. Configure peat layer properties (thickness, hydraulic conductivity, porosity)
4. Prepare climate input data (rainfall, PET)
5. Set model parameters and runtime configuration

### 2. Model Execution
1. Compile source code using provided makefiles
2. Place input files in correct directories
3. Run executable from model_runs directory
4. Monitor simulation progress and output

### 3. Results Processing
1. Extract water table output data
2. Analyse spatial and temporal patterns
3. Create visualisations and summary statistics
4. Generate reports and documentation

## Directory Structure

### `src/`
Contains the DigiBog_Hydro source code and compilation files:
- `DigiBog_Hydro_main.f90` - Main program file
- `DigiBog_Hydro_procs.f90` - Procedures and functions module
- `README.md` - Source code documentation
- `LICENSE.txt` - GNU GPL v3 licence
- `Make_hydro_run`, `Make_hydro_debug` - Makefiles for compilation
- Compiled files (`.o`, `.mod`) and executables

### `parameter_files/`
Model parameter configuration:
- `hydro_parameters.txt` - Scalar parameters for model configuration
  (12 values including grid dimensions, time steps, output intervals)

### `inputs_status_files/`
Field configuration and peat properties:
- `hydro_column_status.txt` - Cell activation status ("on", "off", "diri", "neu")
- `hydro_no_layers.txt` - Number of peat layers per column
- `hydro_thickness.txt` - Layer thickness values (cm)
- `hydro_k.txt` - Hydraulic conductivity values
- `hydro_s.txt` - Drainable porosity values
- `hydro_wt_bc_input.txt` - Initial water table heights/boundary conditions
- `hydro_baltitude.txt` - Base altitude for each column

### `climate_inputs/`
Climate forcing data:
- `hydro_rainfall.txt` - Daily rainfall data
- `hydro_pet.txt` - Daily potential evapotranspiration data

### `model_runs/`
Simulation execution files and run-specific configurations

### `pre_process/`
Scripts and files for model setup and input preparation

### `post_process/`
Scripts and files for results analysis and visualisation

### `documents/`
Project documentation, including this README and implementation guides

## Key Scripts and Their Functions

### Compilation
- **`Make_hydro_run`** - Produces optimised executable for production runs
- **`Make_hydro_debug`** - Produces debug version for troubleshooting
- **Commands:**
  - `make -f Make_hydro_run` - Create optimised executable
  - `make -f Make_hydro_debug` - Create debug executable for use with GDB 
  - `make -f [filename] clean` - Remove executables and object files

### Core Model Functions
- **Water table calculation** - Solves Boussinesq equation
- **AET function** - Calculates actual evapotranspiration based on water table depth
- **Transmissivity calculations** - Determines water movement capacity
- **Boundary condition handling** - Manages Dirichlet and Neumann boundaries
- **Time stepping** - Controls simulation temporal progression

## Running Simulations

### Prerequisites
- GNU Fortran compiler (gfortran)
- GNU Make
- Properly formatted input files

### Input File Requirements

#### Boundary Conditions
- **"on"** - Active simulation cell
- **"off"** - Inactive cell (corner cells are always "off")
- **"diri"** - Dirichlet boundary (constant water height, allows flow)
- **"neu"** - Neumann boundary (no-flow condition)

#### File Formats
- All input files are plain text with numerical data
- Spatial data ordered by x, then y (i.e. all y values for x=1, then x=2, etc.)
- Layer data ordered by x, then y, then z (bottom to top)
- Time series data ordered chronologically
- Scientific notation acceptable (e.g. 1.23e-05)

### Execution Steps
1. **Prepare input files** in respective directories
2. **Compile model:** `make -f Make_hydro_run`
3. **Navigate to model_runs directory**
4. **Run executable:** `./DigiBog_hydro.run` (or debug version)
5. **Monitor output files** for completion

### Output Files
- `hydro_wt_output.txt` - Water table positions relative to base
- `hydro_x_indices.txt` - X coordinates for spatial referencing
- `hydro_y_indices.txt` - Y coordinates for spatial referencing

## Post Processing

### Data Analysis
- **Water table patterns** - Analyse spatial and temporal variations
- **Boundary effects** - Examine influence of different boundary conditions
- **Management scenarios** - Compare pre/post intervention hydrology
- **Statistical summaries** - Calculate means, ranges, trends

### Visualisation
- **Spatial plots** - Water table contour maps
- **Time series** - Temporal changes at specific locations
- **Cross-sections** - Vertical water table profiles
- **Difference maps** - Compare scenario results

### Output Interpretation
- Water table values are relative to impermeable base altitude
- Positive values indicate water table above base
- Negative values indicate water table below base (unusual in peatlands)
- Units are in centimetres above base altitude

### Quality Control
- Check for mass balance closure
- Verify boundary condition behaviour
- Validate against field observations where available
- Assess numerical stability and convergence

## Additional Notes

### Model Limitations
- Water movement only - no peat dynamics
- Assumes homogeneous peat properties within each layer
- Simplified evapotranspiration representation
- No surface water routing

### Best Practices
- Start with simple scenarios before complex setups
- Validate model setup with steady-state tests
- Use appropriate time steps for stability
- Document all assumptions and parameter choices
- Archive input files with results for reproducibility

### Troubleshooting
- Use debug version for identifying numerical issues
- Check input file formats and ordering
- Verify boundary condition setup
- Monitor for unrealistic water table fluctuations
- Consult source code comments for detailed algorithm descriptions

For technical issues or model development questions, contact the code authors listed in the code base.
