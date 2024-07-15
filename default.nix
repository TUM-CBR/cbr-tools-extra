{
  nixpkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz") { },
  netogallo-pkgs ? import (
    fetchGit {
      url = "https://github.com/netogallo/nix-packages.git";
      rev = "bb968081d6cd59189a5a8899fd7d195dbff44bbc";
    }) { inherit nixpkgs; },
}:
let
  python = nixpkgs.python3;
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
  dependencies = import ./requirements.nix { python-pkgs = python.pkgs; python-pkgs-ng = netogallo-pypkgs; };
}

