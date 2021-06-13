with import <nixpkgs> {};

let
  pystan = python38.pkgs.buildPythonPackage rec {
    pname = "pystan";
    version = "2.19.1.1";

    src = python38.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "0f5hbv9dhsx3b5yn5kpq5pwi1kxzmg4mdbrndyz2p8hdpj6sv2zs";
    };

    doCheck = false;

    buildInputs = with python38.pkgs; [
      numpy
      cython
      clikit
    ];
  };

  pythonEnv = python38.withPackages (ps: [
    ps.cython
    ps.numpy
    ps.pandas
    ps.jupyter
    ps.matplotlib
    ps.simpy
  ]);

in mkShell {
  buildInputs = [
    pythonEnv
    pystan
  ];
}


