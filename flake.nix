{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
    flake-parts.url = "github:hercules-ci/flake-parts";
    netogallo-pkgs.url = "github:netogallo/nix-packages?ref=f4bbccc569c22f4c4c186d4ba2f71b31f6f542a7";
  };
  outputs = inputs@{ self, nixpkgs, poetry2nix, flake-parts, netogallo-pkgs, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
	"x86_64-linux"
	"x86_64-darwin"
      ];
      perSystem = { self', pkgs, system, ... }:
      let
	cbr-tools-extra = import ./default.nix {
	  inherit pkgs;
	  netogallo-pkgs = netogallo-pkgs.outputs.legacyPackages.${system};
	}; 
      in
      {
	packages = {
	  inherit cbr-tools-extra;
	  default = cbr-tools-extra;
	};
      };
  };
}
