{pkgs}: let
  # Fetch the source from GitHub
  src = pkgs.fetchFromGitHub {
    owner = "arosl";
    repo = "deco";
    rev = "f2071e8";
    sha256 = "RNpx6VKvBx+/AHL/QFCwvH3+6g4QV4HhhK2d5It54LA=";
  };
in
  with pkgs;
    stdenv.mkDerivation rec {
      pname = "deco";
      version = "1.0.1";

      inherit src;

      installPhase = ''
        mkdir -p $out/bin
        cp deco $out/bin/
      '';

      meta = with lib; {
        description = "Deco - command line tool to create deco profiles";
        license = licenses.bsd3;
        platforms = platforms.unix;
      };
    }
