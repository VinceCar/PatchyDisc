# PatchyDisc
This is a code to perform NVT simulations of patchy particles. 
Several potentials are available.
The code was built based on the [VMMC library](https://github.com/lohedges/vmmc) of Lester Hedges.

## Requirements
- cmake >= 3.11
- [nlohmann json library](https://github.com/nlohmann/json) 

## Compiling PatchyDisc
```
git clone git@github.com:mosayebi/PatchyDisc.git
cd PatchyDisc
mkdir build
cd build
cmake ..
make
```
A successfull compilation will create the executable file `build/PatchyDisc`.

## Running examples

```
cd examples
../build/PatchyDisc input_mixture.json
```
This will run a simulation using the parameters configureed in the `input_mixture.json`.


## IMPORTANT

due to update to the work, examples are currently unavailable. they will be added with the next commit. for the moment refer to the example input below.

## Input file

Below is the contents of a sample input file for simulating a mixture of tri-valent and bi-valent patchy discs.

{
    "initialisation" :
    {
        "mode": "from_random_conf",
        "from_random_conf":
        {
            "types": 2,
            "particle_numbers" : [{"type": 0, "N":5}, {"type": 1, "N": 5}],
            "box": [5.0, 5.0]
        },
        "from_init_conf":
        {
            "types": 2,
            "particle_numbers" : [{"type": 0, "N":1}, {"type": 1, "N":1}],
            "box": [5.0, 5.0],
            "init_conf": "last_conf.xyz",
            "init_patch_conf":"last_states.xvg",
            "restart_step_counter": true
        }
    },


    "topology" :
    {
        "patches": [
            {"type": 0, "nPatches": 2, "angles" : [0.0, 120], "nStates": 1, "probability": [0.0]},
            {"type": 1, "nPatches": 3, "angles" : [0.0, 120.0, 240.0], "nStates": 1, "probability": [0.0]}
        ],
        "pair_coeff" : [
            {"type1":0, "type2":0, "epsilon":8.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.038, "shift": 1.0},
            {"type1":0, "type2":1, "epsilon":8.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.038, "shift": 1.0},
            {"type1":1, "type2":0, "epsilon":8.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.038, "shift": 1.0},
            {"type1":1, "type2":1, "epsilon":8.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.038, "shift": 1.0}
        ]
    },

    "simulation parameters":
    {
        "seed": 844,
        "interaction": "KFOpenSurfPatchyDisc",
        "trajectory": "trajectory.xyz",
        "last_conf": "last_conf.xyz",
        "log_file" : "log.dat",
        "output_every": 1,
        "sweeps": 1e3
    }
}

```
Currently, the following interactions are implemented

- `GaussianPatchyDisc` -- The model is described in this [reference](https://arxiv.org/abs/1607.06626). with the addition to consider if a patch is in an active or inactive state
- `KFSurfOpenPatchyDisc` -- The model is a Kern-Frenkel potential that consider if a patch is in an active or inactive state

----- Not Tested yet
- `GaussianPatchyDiscHR` -- Same as the `GaussianPatchyDisc`, but with a hard core repulsion when particles overlap (i.e. when $r_{ij} \lt \sigma_{ij}$).
- `GaussianPatchyDiscHRSW` -- Same as the `GaussianPatchyDiscHR`, but with a Square-Well poential when particles are within the interaction range (i.e. when $r_{ij} \le {r_{cut}}_{ij}$).

## Important

please note that there is a surface interaction term that is active when particles interact. to turn it off it is necessary to modify in src/gaussian_patchy.cpp the term "surface_interaction"

## Visualisation
The trajectory is saved in XYZ format. You can vissualise it in the visualaisation package of your choice. However, I recommend using [ovito](https://www.ovito.org/).

When importing the trajectory in ovito, please choose the follwing order for the column mappings: Particle type, Position X, Position Y, Orientation Z, Orientation W. You can then load particle shape STL files saved in the `UTILS/particle_shapes/` if you are simulating a standard patchy disc, or make your own particle shape using the python script provided.

To load a particle shape in ovito under the "Data source" panel click on "Particle Types" and choose a particle type. The "Load shape..."  button can then be found in the bottom of the panel.


![Binary mixture of patchy discs](mixture.gif)

## Update
Now it is possible to specify the surface interaction as parameter in the input files. Note that this update breaks the examples listed above, a quick one will be provided in the next month.
