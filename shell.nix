{
  fetchzip ? (import nixpkgs {}).fetchzip,
  nixpkgs ? import (fetchzip {
      url = "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz";
      hash = "sha256-vboIEwIQojofItm2xGCdZCzW96U85l9nDW3ifMuAIdM="; }
    ) {},
  netogallo-pkgs ? import (
    fetchGit {
      url = "https://github.com/netogallo/nix-packages.git";
      rev = "bb968081d6cd59189a5a8899fd7d195dbff44bbc";
    }) { inherit nixpkgs; },
}:
let
  python = nixpkgs.python3;
  netogallo-pypkgs = netogallo-pkgs.python-packages; 
  python-dev = python.withPackages (p:
    import ./requirements.nix { python-pkgs = p; python-pkgs-ng = netogallo-pkgs.python-packages; }
    ++ [ p.pytest p.ipython p.python-lsp-server p.jedi p.pyflakes p.jedi-language-server p.debugpy ]
  );
in
nixpkgs.mkShell {
  name = "cbr-tools-extra";
  packages = with nixpkgs; [
    python-dev
    pyright
    openapi-generator-cli
  ];
}
