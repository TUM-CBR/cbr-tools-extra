{
  python-pkgs,
  python-pkgs-ng
}:
let
  biopython-cbrtools = python-pkgs.biopython.overridePythonAttrs {
    src = fetchGit {
      url = "https://github.com/netogallo/biopython.git";
      # rev = "b609de35cdfef8d437494883f9bb175c63ee58b2"; # 1.8.4
      rev = "a1b846c5fb392d1fdda4cabb5b652264fcd93f95"; # 1.8.3
    };
    doCheck = false;
    # patches = [];
  };
in
with python-pkgs; with python-pkgs-ng; [
    docopt
    igraph
    openpyxl
    open3d-cpu
    pandas
    primer3
    pydantic
    requests
    sqlalchemy
    biopython-cbrtools
]
