let
  nixpkgs = import <nixpkgs> {};
  inherit (nixpkgs) callPackage;
  inherit (nixpkgs.python3.pkgs) buildPythonPackage;
  embree3 = callPackage ./embree/embree3.nix { };
in
{
  inherit embree3;
  open3d-cpu = callPackage ./open3d/open3d-cpu.nix { inherit buildPythonPackage embree3; };
}

