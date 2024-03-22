#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Check operating system and download appropriate installation file
case "`uname -s`" in
    Linux*)
        echo "Recognized Linux as OS (assuming x86_64), installing..."
        wget https://github.com/Kitware/CMake/releases/download/v3.19.1/cmake-3.19.1-Linux-x86_64.tar.gz -O cmake.tar.gz
        # Untar the file
        tar -zxvf cmake.tar.gz
        # Rename directory
        mv cmake-3.19.1-Linux-x86_64 cmake
        # Remove the tar ball
        rm cmake.tar.gz
    ;;
    Darwin*)
        echo "Recognized macOS as OS (assuming x86_64), installing..."
        curl -o cmake.tar.gz https://github.com/Kitware/CMake/releases/download/v3.19.1/cmake-3.19.1-Darwin-x86_64.tar.gz
        # Untar the file
        tar -zxvf cmake.tar.gz
        # Rename directory
        mv cmake-3.19.1-Darwin-x86_64/CMake.app/Contents cmake
        # Remove the tar ball
        rm cmake.tar.gz
    ;;
    *)
        echo "Could not recognize OS. Aborting."
        return
esac

# Get absolute path to binary folder
thisPath=`pwd`"/cmake/bin/"

# Prompt user with the relevant path
echo ""
echo "Don't forget to add the following path to your shell's PATH variable:"
echo $thisPath 
echo ""
echo "This can be done e.g. with:"
echo "export PATH=\""$thisPath"\":\${PATH}"
