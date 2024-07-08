{
  description = ''
  This flake allows running cbr-tools-extra. It is a program used to support cbr-tools by bundling
  various bioniformatics programs and implementing various bioinformatics algorithms.
  '';

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
    flake-parts.url = "github:hercules-ci/flake-parts";
    netogallo-pkgs = {
      url = "github:netogallo/nix-packages?ref=bb968081d6cd59189a5a8899fd7d195dbff44bbc";
      inputs.nixpkgs.follows = "nixpkgs";
    };
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
	  nixpkgs = pkgs;
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
