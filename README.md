# cbr-tools-extra
Additional components for CBR tools which rely on libraries not usually shipped with a PyMol distribution.

## Setup

### Linux

The easiest way to setup this program for use on Linux (with cbr-tools) is using a virtual environment. Breif description is provided below:

1. Ensure your terimal is scoped to the directory containing this README.md file.
2. Ensure you have [Python](https://www.python.org/) and [virtualenv](https://virtualenv.pypa.io/en/latest/user_guide.html) installed.
3. Create a virtualenv. You must name it _.virtualenv as that is the name used by the _cbrtools_ script. Activate it afterwards.
```
virtualenv .virtualenv
source .virtualenv/bin/activate
```
4. Install the necessary dependencies
```
pip install -r requirements.base.txt
```
5. Make the _cbrtools_ script executable
```
chmod +x cbrtools
```
6. The program is ready to run, you can test it with:
```
./cbrtools --help
```
7. If you are using this program in conjunction with [CBR Bioinformatics Tools](https://github.com/TUM-CBR/pymol-plugins), the _cbrtools_ script must be in your path. The easiest way to achive this is to add this directory to your path by editing your _bashrc_ file. A one liner to do that is:
```
echo "export PATH=`pwd`:\$PATH" >> ~/.bashrc
```

### Windows

 For Windows users we advise that you download the latest [nightly build](https://github.com/TUM-CBR/cbr-tools-extra/releases). If you are using this in conjunction with [CBR Bioinformatics Tools](https://github.com/TUM-CBR/pymol-plugins), no setup should be necessary as this tools are bundled in the [nightly builds of CBR Bioinformatics Tools](https://github.com/TUM-CBR/pymol-plugins/releases).