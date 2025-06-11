# ESIGMAPy: a Python package to generate `ESIGMAHM` waveforms

`ESIGMAHM` is an eccentric, aligned-spin, inspiral-merger-ringdown (IMR) waveform model with higher-order modes. It is composed of two pieces:

* **Inspiral piece (called `InspiralESIGMAHM`):** The inspiral piece comes from a combination of post-Newtonian theory, self-force, and black hole perturbation theory. It is currently implemented in a [private fork](https://git.ligo.org/kaushik.paul/lalsuite/-/tree/enigma_spins_v2023?ref_type=heads) of `LALSuite` (interested users are welcome to write to the developers for access at esigmahm@icts.res.in).  
* **Plunge-merger-ringdown piece**: Assuming moderate starting eccentricities that will decay by the late inspiral, we use the quasi-circular NR surrogate `NRSur7dq4` for the plunge-merger-ringdown piece for `ESIGMAHM`. This requires the `NRSur7dq4`'s [data file](https://git.ligo.org/lscsoft/lalsuite-extra/-/blob/master/data/lalsimulation/NRSur7dq4.h5) to be downloaded. 

  (_Note:_ We also allow other `LALSuite` waveform models to be used as the plunge-merger-ringdown piece; see the optional argument `merger_ringdown_approximant` in the generation functions `get_imr_esigma_waveform` and `get_imr_esigma_modes` [here](https://github.com/gwnrtools/esigmapy/blob/master/esigmapy/generator.py). However, the default and the most tested choice is `NRSur7dq4`)

The full IMR waveform `ESIGMAHM` is then produced by smoothly attaching the inspiral piece `InspiralESIGMAHM` to the plunge-merger-ringdown piece `NRSur7dq4`. This attachment is done via `ESIGMAPy`.

## Installing `InspiralESIGMAHM`

* **Getting the source code:** The `LALSuite` fork containing the implementation of `InspiralESIGMAHM` is currently private, but interested users are welcome to write to the developers for access at esigmahm@icts.res.in. 

  Clone this [`LALSuite` fork](https://git.ligo.org/kaushik.paul/lalsuite/-/tree/enigma_spins_v2023?ref_type=heads) and checkout the relevant commit identified by the tag `ESIGMAHMv1` which contains the current stable implementation of the inspiral piece of `ESIGMAHM`.

  ```
  git clone https://git.ligo.org/kaushik.paul/lalsuite.git
  cd lalsuite
  git checkout ESIGMAHMv1
  ``` 
* **Installing the code:**
  - Activate your `conda` environment. Make sure that the `swig` version in this environment is below `4.2.1` (you can check this by running `conda list swig`). If not, install its version `4.2.0` by running `conda install -c conda-forge swig=4.2.0`. 
  - Now choose/create a directory where you want to install ESIGMA. Let the absolute path of this directory be `/path/to/esigmahm`.
  - Go back inside the above cloned `LALSuite` fork, and sequentially run the following commands: 
    
    ```
    ./00boot
    ./configure --prefix="/path/to/esigmahm" --enable-swig-python --disable-laldetchar --disable-lalpulsar --disable-lalapps --enable-mpi=yes --enable-hdf5 CFLAGS="-Wno-error" CXXFLAGS="-Wno-error" CPPFLAGS="-Wno-error" 
    make
    make install
    ```
* **Configuring `conda` to source this `LALSuite` fork automatically on activating the `conda` environment**
  - With your `conda` environment activated, run the following commands:
    ```
    cd $CONDA_PREFIX/etc/conda/activate.d
    ln -s /path/to/esigmahm/etc/lalsuite-user-env.sh
    ```

## Installing `NRSur7dq4`
* Download the `NRSur7dq4` [data file](https://git.ligo.org/lscsoft/lalsuite-extra/-/blob/master/data/lalsimulation/NRSur7dq4.h5) in some directory. Let's say the absolute path to this directory is `/path/to/NRSur7dq4`.
* Append the path of this directory to the shell environment variable `LAL_DATA_PATH` by running: `export LAL_DATA_PATH="$LAL_DATA_PATH:/path/to/NRSur7dq4"`
* To avoid performing the above step in every new terminal session, either add the above command to your `.bashrc` file, or follow the instructions [here](http://gitlab.icts.res.in/akash.maurya/Installation-instructions/wikis/conda-tricks), replacing `PYTHONPATH` with `LAL_DATA_PATH`, to set this environment variable automatically on activating your `conda` environment.

## Installing `ESIGMAPy`
* Activate your `conda` environment and install `ESIGMAPy` by running: `pip install esigmapy`.

***
## Trying out `ESIGMAHM`
If everything goes fine, you should be able to generate `ESIGMAHM` waveforms. The instructions to do so and the various functionalities that `ESIGMAHM` offers are detailed in [this tutorial notebook](https://github.com/gwnrtools/esigmapy/blob/master/notebooks/ESIGMA_generation.ipynb). 

***
## Citation
If you use `ESIGMAHM` in your work, please consider citing it: 

Paul et. al., _"ESIGMAHM: An Eccentric, Spinning inspiral-merger-ringdown waveform model with Higher Modes for the detection and characterization of binary black holes"_, DOI:[https://doi.org/10.1103/PhysRevD.111.084074](https://doi.org/10.1103/PhysRevD.111.084074), arXiv:[2409.13866](https://arxiv.org/abs/2409.13866) (2024)
```
@article{Paul:2024ujx,
    author = "Paul, Kaushik and Maurya, Akash and Henry, Quentin and Sharma, Kartikey and Satheesh, Pranav and Divyajyoti and Kumar, Prayush and Mishra, Chandra Kant",
    title = "{Eccentric, spinning, inspiral-merger-ringdown waveform model with higher modes for the detection and characterization of binary black holes}",
    eprint = "2409.13866",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.111.084074",
    journal = "Phys. Rev. D",
    volume = "111",
    number = "8",
    pages = "084074",
    year = "2025"
}
```
`ESIGMAHM` is built on the `ENIGMA` framework, which was developed in arXiv:[1609.05933](https://arxiv.org/abs/1609.05933), arXiv:[1711.06276](https://arxiv.org/abs/1711.06276), arXiv:[2008.03313](https://arxiv.org/abs/2008.03313). Thus, in addition to citing `ESIGMAHM`, please consider citing these works related to ENIGMA as well. 

***
## 📬 Contact Us  
If you have any questions, issues, or suggestions regarding the model, feel free to reach out to us at esigmahm@icts.res.in!
