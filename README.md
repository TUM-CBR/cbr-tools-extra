# CBR Tools Extras

CBR Tools Extras is a binary package that contains supporting functionality for [CBR Tools](https://github.com/TUM-CBR/pymol-plugins). We advise you look at that project first if you are not familiar with it. Even though this project is meant to be used in conjunction with [CBR Tools](https://github.com/TUM-CBR/pymol-plugins), it can also work as a standalone command-line application.


## Installation

You can get a packaged binary from the [Releases](https://github.com/TUM-CBR/cbr-tools-extra/releases) section. If you instead wish to create the binary yourself, the process is described below.

1. **System Dependencies**: This program needs [Python](https://www.python.org/) to work. On Linux, you will also need [tcl/tk](https://tcl.tk/doc/).
2. **Clone the Repository**: Use git to get a copy of this repository: `git clone git@github.com:TUM-CBR/cbr-tools-extra.git; cd cbr-tools-extra/`
3. **Python Dependencies**: Use pip to install the dependencies of this project: `pip install -r requirements.base.txt`
4. **Project Ready**: Try running the project with `python cbrtools.py --help`
5. **Testing**: You can run this project's tests with [pytest](https://docs.pytest.org/en/8.0.x/contents.html) by simply running pytest in the root directory: `pytest`
6. **Create a package**: You can use [PyInstaller](https://pyinstaller.org/en/stable/) to package this project. A conveniece script is provided at _scripts/bundle.sh_ which generates a binary inside the _dist_ directory: `sh scripts/bundle.sh`

# Usage and Features

This is a command line program meant to be used with [CBR Tools](https://github.com/TUM-CBR/). To summarize some of the moduels:

1. **cascades**: You can look at this module by running `python cbrtools.py cascades --help`. This module provides the underlying functionality for [Cascade BLAST](https://github.com/TUM-CBR/pymol-plugins/wiki/Cascade-BLAST)
2. **cavities**: You can look at this module by running `python cbrtools.py cavities --help`. This module provides the underlying functionality for [Prot o' Dentist](https://github.com/TUM-CBR/pymol-plugins/wiki/Prot-o'-Dentist)

## License

CBR Tools Extras is released under the "GNU General Public License Version 2". This license allows you to freely use and re-distribtue the software (even for commercial and for-profit reasons) with the condition that any modifications or improvements are made available to the general public under the same license.
