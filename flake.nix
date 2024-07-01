{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
    poetry2nix = {
      url = "github:nix-community/poetry2nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
  outputs = inputs@{ self, nixpkgs, poetry2nix, flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
	"x86_64-linux"
	"x86_64-darwin"
      ];
      perSystem = { self', pkgs, system, ... }:
      let
	poetry = poetry2nix.lib.mkPoetry2Nix { inherit pkgs; };
	inherit (poetry) mkPoetryApplication;
	withPythonLibs =
	  { prev, final, pkg, libs, ... }: prev.${pkg}.overridePythonAttrs(
	      old: {
		buildInputs =
		  (old.buildInputs or [])
		  ++ map (lib: prev.${lib}) libs;
	      }
	  );
	python = pkgs.python311;
      in
      {
	packages = {
	  cbrtools = mkPoetryApplication {
	    inherit python;
	    projectDir = self;
	    overrides = poetry.overrides.withDefaults (final: prev: {
	      primer3-py = withPythonLibs { inherit prev final; pkg = "primer3-py"; libs = [ "setuptools" ]; };
	      pyquaternion = prev.pyquaternion.override {
		preferWheel = true;
	      };
	      open3d-cpu = python.pkgs.numpy;
	    });
	  };
	  default = self'.packages.cbrtools;
	};
	devShells = {
	  default = pkgs.mkShell {
	    inputsFrom = [ self'.packages.cbrtools ];
	  };
	  poetry = pkgs.mkShell {
	    packages = [ pkgs.poetry.override { python3 = pkgs.python312; } ];
	  };
	};
      };
    };
}
