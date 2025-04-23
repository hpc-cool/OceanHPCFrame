# OceanHPCFrame
# OceanHPCFrame

OceanHPCFrame is a parallel framework for accelerating ocean simulation models. It includes four progressively optimized versions to evaluate the performance impact of vectorization, multithreading, and multi-node parallelism.

## Build Instructions

Use the `make` command to compile each version:

- `make v1`: Baseline implementation 
- `make v2`: Vectorized version
- `make v3`: Multi-threaded version
- `make v4`: Fully optimized version

## Run Instructions

Use the corresponding script to run:

- `./run_v1`:
- `./run_v2`: 
- `./run_v3`:
- `./run_v4 <n>`

## Configuration

- Simulation parameters are configured in `nwamctl.ini`, including grid resolution (`GRDSZX`, `GRDSZY`)
- Input topography data is provided in `etop5.nc`

## Environment Requirements

This project was tested on the LS Pilot System with the following software:

- **MPI** v4.1.6  
- **netCDF** v4.9.2  
- **Clang** v17.0.6  
- **Flang** v17.0.6  
- **Sdma** v1.0  
- **Memkind** v1.0

Please make sure your environment includes these or compatible versions for execution.

---

Feel free to open an issue if you have any questions.