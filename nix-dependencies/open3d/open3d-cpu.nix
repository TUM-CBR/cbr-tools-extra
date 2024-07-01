{
  lib,
  pkgs,
  buildPythonPackage,
  embree3
}:
import ./open3d-common.nix {
  inherit lib pkgs buildPythonPackage embree3;
}
