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
  python-dev = python.withPackages (p:
    import ./requirements.nix { python-pkgs = p; python-pkgs-ng = netogallo-pkgs.python-packages; }
    ++ [ p.pytest p.ipython ]
  );
in
nixpkgs.mkShell {
  name = "cbr-tools-extra";
  packages = [
    python-dev
  ];
}
