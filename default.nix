{
  pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz") { },
  netogallo-pkgs ? import (fetchGit { url = "https://github.com/netogallo/nix-packages.git"; rev = "f4bbccc569c22f4c4c186d4ba2f71b31f6f542a7"; }) { nixpkgs = pkgs; }
}:
let
  python = pkgs.python3;
  netogallo-pypkgs = netogallo-pkgs.python-packages; 
  biopython = python.pkgs.biopython.overridePythonAttrs {
    src = fetchGit { url = "https://github.com/netogallo/biopython.git"; rev = "a95ae4130580d09107882ee6bdbc159d4803b122"; };
  };
in
python.pkgs.buildPythonPackage {
  pname = "cbr-tools-extra";
  version = "0.1.0";
  src = ./.;
  pyproject = true;
  dependencies = with python.pkgs; with netogallo-pypkgs; [
    open3d-cpu
    docopt
    igraph
    openpyxl
    open3d-cpu
    pandas
    primer3
    pydantic
    requests
    sqlalchemy
    biopython
  ];
}

