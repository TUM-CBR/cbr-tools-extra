{ pkgs, lib, ... }:
let
  mkDerivation = pkgs.stdenv.mkDerivation;
  inherit (pkgs.python311Packages) virtualenv;
  inherit (pkgs) git;
in
mkDerivation {
  name = "cbr-tools-extra";
  packages = with pkgs; [ virtualenv ];
  src = ./.;
  depsBuildBuild = [ git virtualenv ];
  buildPhase = ''
    virtualenv $out/python
    $out/python/bin/pip install -r requirements.nix.txt
    echo "echo 'cbrtoosl'" > $out/bin/cbrtools
  '';
}

