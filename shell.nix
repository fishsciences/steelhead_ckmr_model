with import <nixpkgs> {};

let
  tqdm = python39.pkgs.buildPythonPackage rec {
    pname = "tqdm";
    version = "4.66.1";
    format = "wheel";

    src = python39.pkgs.fetchPypi {
      inherit pname version format;
      sha256 = "d302b3c5b53d47bce91fea46679d9c3c6508cf6332229aa1e7d8653723793386";
      python = "py3";
      dist = "py3";
    };

    doCheck = false;
  };

  xarray-einstats = python39.pkgs.buildPythonPackage rec {
    pname = "xarray-einstats";
    version = "0.6.0";
    format = "pyproject";

    src = python39.pkgs.fetchPypi {
      inherit version;
      pname = "xarray_einstats";
      sha256 = "ace90601505cfbe2d374762e674557ed14e1725b024823372f7ef9fd237effad";
    };

    doCheck = false;

    propagatedBuildInputs = with python39.pkgs; [
      numpy
      cython  
      scipy
      xarray
      flit
      flit-core
      packaging
      typing-extensions
      netcdf4
    ];
  };

  stanio = python39.pkgs.buildPythonPackage rec {
    pname = "stanio";
    version = "0.3.0";
    format = "wheel";

    src = python39.pkgs.fetchPypi {
      inherit pname version format;
      sha256 = "2d34b5ebe9ad8fcd137437209bf4b53846d88dbe933441aca5d83fd32f9b0c7e";
      python = "py3";
      dist = "py3";
    };

    doCheck = false;

    propagatedBuildInputs = with python39.pkgs; [
      numpy
    ];
  };

  arviz = python39.pkgs.buildPythonPackage rec {
    pname = "arviz";
    version = "0.16.1";

    src = python39.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "35bab9072f66f5a8204d2a71911d09ce05056c177f1a780de53efa2714c27575";
    };

    doCheck = false;

    propagatedBuildInputs = with python39.pkgs; [
      numpy
      cython  
      pandas
      matplotlib
      xarray
      xarray-einstats
      h5netcdf
    ];
  };

  cmdstanpy = python39.pkgs.buildPythonPackage rec {
    pname = "cmdstanpy";
    version = "1.2.0";

    src = python39.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "bdf55ab76f9eea01763b8990a56ff55d03e69ec31d9613fdbbe4c452126ff1bb";
    };

    doCheck = false;

    propagatedBuildInputs = with python39.pkgs; [
      numpy
      cython  
      ujson
      pandas
      tqdm
      stanio
    ];
  };

  pysimdjson = python39.pkgs.buildPythonPackage rec {
      pname = "pysimdjson";
      version = "5.0.2";

      src = python39.pkgs.fetchPypi {
          inherit pname version;
          sha256 = "83010f07f9ca38e4557b61860acfeb0a897b416f06f73182ffaffa94bdb7394d";
      };

      doCheck = false;

      propagatedBuildInputs = with python39.pkgs; [
        pybind11
      ];
  };

  httpstan = python39.pkgs.buildPythonPackage rec {
    pname = "httpstan";
    version = "4.10.1";
    format = "wheel";

    src = python39.pkgs.fetchPypi {
      inherit pname version format;
      sha256 = "821469883091523fef2c6d806bbfdb43efa07ec0df662ba58ccda9e773224640";
      dist = "cp311";
      python = "cp311";
      abi = "cp311";
      platform = "manylinux_2_17_x86_64.manylinux2014_x86_64";
    };

    doCheck = false;

    propagatedBuildInputs = with python39.pkgs; [
      numpy
      cython
      pandas
      clikit
      aiohttp
      pysimdjson
      marshmallow
      appdirs
      webargs
    ];
  };

  pythonEnv = python39.withPackages (ps: [
    ps.cython
    ps.numpy
    ps.pandas
    ps.openpyxl
    #ps.jupyterlab
    #ps.altair
    ps.matplotlib
    #ps.ipywidgets
  ]);

in mkShell {
  buildInputs = [
    pythonEnv
    cmdstanpy
    arviz
    tqdm
  ];
}


