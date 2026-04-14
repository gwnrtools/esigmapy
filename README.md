# ESIGMAPy: a Python package to generate ESIGMA waveforms

`ESIGMAPy` provides access to two eccentric waveform models:

| Model | Description |
|---|---|
| `ESIGMAHM` | Eccentric, aligned-spin IMR waveform with higher modes (arXiv:[2409.13866](https://arxiv.org/abs/2409.13866)) |
| `ESIGMASur` | Surrogate model of the (2,2) GW-mode of `ESIGMAHM` (arXiv:[2510.00116](https://arxiv.org/abs/2510.00116))|

---

## :rocket: Quick Start

Install via:

```bash
pip install git+https://github.com/gwnrtools/esigmapy.git
```
Then to use:
 * `ESIGMAHM`: requires custom `LALSuite` installation and `NRSur7dq4` surrogate data download ([installation instructions below](https://github.com/gwnrtools/esigmapy/tree/master?tab=readme-ov-file#blue_square-esigmahm)).
 * `ESIGMASur`: requires downloading `ESIGMASur` surrogate data ([installation instructions below](https://github.com/gwnrtools/esigmapy/tree/master?tab=readme-ov-file#green_square-esigmasur))
 
Usage instructions at:
 * `ESIGMAHM` [tutorial notebook](https://github.com/gwnrtools/esigmapy/blob/master/notebooks/ESIGMA_tutorial.ipynb)
 * `ESIGMASur` [tutorial notebook](https://github.com/gwnrtools/esigmapy/blob/master/notebooks/ESIGMASur_tutorial.ipynb)
***
## :blue_square: ESIGMAHM

`ESIGMAHM` is an eccentric, aligned-spin, inspiral-merger-ringdown (IMR) waveform model with higher-order modes. It is composed of two pieces:

* **Inspiral piece (`InspiralESIGMAHM`):** It comes from a combination of post-Newtonian theory, self-force, and black hole perturbation theory.
* **Plunge-merger-ringdown piece**: Assuming moderate starting eccentricities that decay by the late inspiral, we use the quasi-circular (QC) NR surrogate `NRSur7dq4` for the plunge-merger-ringdown piece for `ESIGMAHM`.
  
  (**Note:** We also allow other QC `LALSuite` waveform models to be used as the plunge-merger-ringdown piece; see the optional argument `merger_ringdown_approximant` in the generation functions `get_imr_esigma_waveform` and `get_imr_esigma_modes` [here](https://github.com/gwnrtools/esigmapy/blob/master/esigmapy/generator.py). However, the default and the most tested choice is `NRSur7dq4`.)

The full IMR waveform `ESIGMAHM` is produced by smoothly attaching the inspiral piece `InspiralESIGMAHM` to the plunge-merger-ringdown piece `NRSur7dq4`. This attachment is facilitated by `ESIGMAPy`.

Using `ESIGMAHM` therefore requires installing 1. `InspiralESIGMAHM`, 2. `NRSur7dq4`, and 3. `ESIGMAPy`.

### Installing `InspiralESIGMAHM`

* **Getting the source code:** The `LALSuite` fork containing the implementation of `InspiralESIGMAHM` is currently private, but interested users are welcome to write to the developers for access at esigmahm@icts.res.in. 

  Clone the [`LALSuite` fork](https://git.ligo.org/kaushik.paul/lalsuite/-/tree/enigma_spins_v2023?ref_type=heads) and checkout the commit with tag `ESIGMAHMv1`:

  ```bash
  git clone https://git.ligo.org/kaushik.paul/lalsuite.git
  cd lalsuite
  git checkout ESIGMAHMv1
  ``` 
* **Installation:**
  - It is advised to create a conda environment using [the `igwn` yaml file](https://computing.docs.ligo.org/conda/environments/igwn-py311/) and install the above mentioned fork of `LALSuite` inside it to minimize the dependency issues. Please remove `_x86_64-microarch-level=3=3_haswell` from the yaml file in case you get errors related to microarchitecture mismatch.  
  - Activate your `conda` environment. Make sure that the `swig` version in this environment is below `4.2.1` (you can check this by running `conda list swig`). If not, install its version `4.2.0` by running `conda install -c conda-forge swig=4.2.0`. 
  - Create a directory where you want to install ESIGMA. Let the absolute path of this directory be `/path/to/esigmahm`.
  - Go back inside the above cloned `LALSuite` fork, and run the following commands: 
    
    ```bash
    ./00boot
    ./configure --prefix="/path/to/esigmahm" --enable-swig-python --disable-laldetchar --disable-lalpulsar --disable-lalapps --enable-mpi=yes --enable-hdf5 CFLAGS="-Wno-error" CXXFLAGS="-Wno-error" CPPFLAGS="-Wno-error" 
    make
    make install
    ```
* **Configuring `conda` to source this `LALSuite` fork automatically on activating the `conda` environment**
  - With your `conda` environment activated, run the following commands:
    ```bash
    cd $CONDA_PREFIX/etc/conda/activate.d
    ln -s /path/to/esigmahm/etc/lalsuite-user-env.sh
    ```

### Installing `NRSur7dq4`
* Download the `NRSur7dq4` [data file](https://git.ligo.org/lscsoft/lalsuite-extra/-/blob/master/data/lalsimulation/NRSur7dq4.h5) in some directory (for `LALSuite` versions >= 7.25, download [this data file](https://git.ligo.org/waveforms/software/lalsuite-waveform-data/-/blob/main/waveform_data/NRSur7dq4_v1.0.h5?ref_type=heads) instead). Let's say the absolute path to this directory is `/path/to/NRSur7dq4`.
* Append the path of this directory to the shell environment variable `LAL_DATA_PATH` by running: `export LAL_DATA_PATH="$LAL_DATA_PATH:/path/to/NRSur7dq4"`
* To avoid performing the above step in every new terminal session, either add the above command to your `.bashrc` file, or follow the instructions [here](http://gitlab.icts.res.in/akash.maurya/Installation-instructions/wikis/conda-tricks), replacing `PYTHONPATH` with `LAL_DATA_PATH`, to set this environment variable automatically on activating your `conda` environment.

### Installing `ESIGMAPy`
* Activate your `conda` environment and install `ESIGMAPy` by running: `pip install esigmapy`.

***
### Trying out `ESIGMAHM`
The instructions to generate `ESIGMAHM` waveforms and the various functionalities that it offers are detailed in [this tutorial notebook](https://github.com/gwnrtools/esigmapy/blob/master/notebooks/ESIGMA_tutorial.ipynb). 

### Trying out `ESIGMAHM` via `gwsignal` 
`ESIGMAHM` can now also be accessed via [`gwsignal`](https://waveforms.docs.ligo.org/reviews/lalsuite/lalsimulation/gwsignal/index.html). See [this notebook](https://github.com/gwnrtools/esigmapy/blob/cfe881a6052f1e9a7c3e066d9fd32462a3af2c89/notebooks/esigma_gwsignal.ipynb) for usage instructions.

***

## :green_square: ESIGMASur

`ESIGMASur` is a fast time-domain surrogate model of the (2,2)-mode of `ESIGMAHM`.

### Installation
* Install `ESIGMAPy` by running:
  ```bash
  pip install git+https://github.com/gwnrtools/esigmapy.git
  ```
  It **does not** require the custom LALSuite installation above. 

  **Note:** `ESIGMASur` currently supports Python <= 3.12, owing to a compatibility restriction in [`tpi-splines`](https://github.com/mpuerrer/TPI), a dependency of `ESIGMASur`.
    
* Download the surrogate data files of `ESIGMASur`, which can be found [on the repo](https://github.com/gwnrtools/esigmapy/tree/master/esigmapy/surrogate/data). Next, set the shell environment variable `ESIGMASUR_DATA_PATH` to the directory where you keep these surrogate data files by running: `export ESIGMASUR_DATA_PATH="/path/to/ESIGMASur"`.

* Using the inspiral-only surrogate `InspiralESIGMASur` would not require any further dependencies. However, its hybridized IMR version `IMRESIGMASur` will require downloading the NR surrogate data file ([installation instructions above](https://github.com/gwnrtools/esigmapy/tree/master#installing-nrsur7dq4)). 

### Trying out `ESIGMASur`
The usage instructions and the various functionalities of `ESIGMASur` are detailed in [this tutorial notebook](https://github.com/gwnrtools/esigmapy/blob/master/notebooks/ESIGMASur_tutorial.ipynb).

***
## 📖 Citation
* If you use `ESIGMAHM` in your work, please cite:
  
  Paul et al., _"ESIGMAHM: An Eccentric, Spinning inspiral-merger-ringdown waveform model with Higher Modes for the detection and characterization of binary black holes"_, [PhysRevD.111.084074](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.111.084074), arXiv:[2409.13866](https://arxiv.org/abs/2409.13866) (2024)

  ```bibtex
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
  `ESIGMAHM` is built on the `ENIGMA` framework, which was developed in arXiv:[1609.05933](https://arxiv.org/abs/1609.05933), arXiv:[1711.06276](https://arxiv.org/abs/1711.06276), arXiv:[2008.03313](https://arxiv.org/abs/2008.03313). In addition to citing `ESIGMAHM`, please consider citing them as well. 

* If you use `ESIGMASur` in your work, please cite:

  Maurya et al., _"Chase Orbits, not Time: A Scalable Paradigm for Long-Duration Eccentric Gravitational-Wave Surrogates"_, arXiv:[2510.00116](https://arxiv.org/abs/2510.00116) (2025)

  ```bibtex
  @article{Maurya:2025shc,
      author = "Maurya, Akash and Kumar, Prayush and Field, Scott E. and Mishra, Chandra Kant and Nee, Peter James and Paul, Kaushik and Pfeiffer, Harald P. and Ravichandran, Adhrit and Varma, Vijay",
      title = "{Chase Orbits, not Time: A Scalable Paradigm for Long-Duration Eccentric Gravitational-Wave Surrogates}",
      eprint = "2510.00116",
      archivePrefix = "arXiv",
      primaryClass = "gr-qc",
      journal = {},
      month = "9",
      year = "2025",
  }
  ```

***
## :mailbox_with_mail: Contact Us  
If you have any questions, issues, or suggestions regarding the model, feel free to reach out to us at esigmahm@icts.res.in!
