#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Check operating system and download appropriate installation file
case "`uname -s`" in
    Linux*)
        echo "Recognized Linux as OS (assuming x86_64), installing..."
        wget https://gmsh.info/bin/Linux/gmsh-4.7.1-Linux64-sdk.tgz -O gmsh.tgz
        # Untar the file
        tar -zxvf gmsh.tgz
        # Rename directory
        mv gmsh-4.7.1-Linux64-sdk gmsh
        # Remove the tar ball
        rm gmsh.tgz
    ;;
    Darwin*)
        echo "Recognized macOS as OS (assuming x86_64), installing..."
        wget https://gmsh.info/bin/MacOSX/gmsh-4.7.1-MacOSX-sdk.tgz -O gmsh.tgz
        # Untar the file
        tar -zxvf gmsh.tgz
        # Rename directory
        mv gmsh-4.7.1-MacOSX-sdk gmsh
        # Remove the tar ball
        rm gmsh.tgz
    ;;
    *)
        echo "Could not recognize OS. Aborting."
        return
esac
