with import <nixpkgs> {};

let
  xarray-einstats = python38.pkgs.buildPythonPackage rec {
    pname = "xarray-einstats";
    version = "0.3.0";
    format = "pyproject";

    src = python38.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "13h9fbbz43y9diczv2b7xm4y3dv38zrnj5lz8qkr6iqqa8a7q8c1";
    };

    doCheck = false;

    propagatedBuildInputs = with python38.pkgs; [
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

  arviz = python38.pkgs.buildPythonPackage rec {
    pname = "arviz";
    version = "0.12.1";

    src = python38.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "06ghp01vxnplsl3a9zy5fyqc69r2n2k9cg5nda7g228rqnn0xn2p";
    };

    doCheck = false;

    propagatedBuildInputs = with python38.pkgs; [
      numpy
      cython  
      pandas
      matplotlib
      xarray
      xarray-einstats
    ];
  };

  cmdstanpy = python38.pkgs.buildPythonPackage rec {
    pname = "cmdstanpy";
    version = "0.9.76";

    src = python38.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "17y8mkqgxfcbc7l2fcyx4nqimrx2jzymyp6p74rfywfq9d22zk2n";
    };

    doCheck = false;

    propagatedBuildInputs = with python38.pkgs; [
      numpy
      cython  
      ujson
      pandas
    ];
  };

  pysimdjson = python38.pkgs.buildPythonPackage rec {
      pname = "pysimdjson";
      version = "3.2.0";

      src = python38.pkgs.fetchPypi {
          inherit pname version;
          sha256 = "0xqlxjr2s369i7pl3gv7ryba2k5zabsin2dw3mv6f8vm844slfv4";
      };

      doCheck = false;

      propagatedBuildInputs = with python38.pkgs; [
        pybind11
      ];
  };

  pystan = python38.pkgs.buildPythonPackage rec {
    pname = "pystan";
    version = "3.5.0";

    src = python38.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "0aa9a4szwi0cch5p4jhznv24m7wpl2a4jv90b6pw1dx5f7873187";
    };

    doCheck = false;

    propagatedBuildInputs = with python38.pkgs; [
      numpy
      cython
      pandas
      clikit
      aiohttp
      pysimdjson
      httpstan
    ];
  };

  httpstan = python38.pkgs.buildPythonPackage rec {
    pname = "httpstan";
    version = "4.8.1";
    format = "wheel";

    src = python38.pkgs.fetchPypi {
      inherit pname version format;
      sha256 = "1lv5gdiafi4w07kgsbnmfcp1fn2imdax2h6q8ncb732zih52yqww";
      dist = "cp38";
      python = "cp38";
      abi = "cp38";
      platform = "manylinux_2_17_x86_64.manylinux2014_x86_64";
    };

    doCheck = false;

    propagatedBuildInputs = with python38.pkgs; [
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

  pythonEnv = python38.withPackages (ps: [
    ps.cython
    ps.numpy
    ps.pandas
    ps.openpyxl
    ps.jupyterlab
    ps.altair
    ps.matplotlib
    ps.tqdm
    ps.ipywidgets
  ]);

in mkShell {
  buildInputs = [
    pythonEnv
    pystan
    cmdstanpy
    arviz
  ];
}


