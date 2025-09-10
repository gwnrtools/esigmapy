# pycbc_esigma_plugin
PyCBC plugin to call ESIGMA waveforms

This plugin can be installed after installation of usual PyCBC in order to call ESIGMA waveforms in PyCBC framework. 
Note: This requires gwnrtools to be installed.

## Installation Instructions
### Install gwnr package
```
git clone https://github.com/gwnrtools/gwnrtools.git
cd gwnr
pip install -r requirements.txt
pip install .
```

### Install pycbc_esigma_plugin
Install the plugin
```
git clone https://github.com/divyajyoti09/pycbc_esigma_plugin.git
cd pycbc_esigma_plugin
pip install .
```
Set `$LAL_DATA_PATH` which contains essential waveforms files such as `NRSur7dq4.h5`. This is required to call IMRESIGMA waveforms
```
export LAL_DATA_PATH=/path/to/lal_data_files
```
Run the post-install-script. This script modifies `gwnr/waveform/esigma_utils.py` to include the path to `$LAL_DATA_PATH` specified above
```
pycbc-esigma-post-install
```
