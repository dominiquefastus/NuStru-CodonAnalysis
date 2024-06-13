#!/bin/bash

# Function to check if a program is installed
is_installed() {
    if command -v $1 >/dev/null 2>&1; then
        echo "$1 already installed: $(which $1)"
    else
        echo "$1 not installed"
        return 1
    fi
}

# Function to install a program
install_program() {
    if is_installed $1; then
        return
    fi

    if [[ "$unamestr" == 'Linux' ]]; then
        sudo apt install $1
    elif [[ "$unamestr" == 'Darwin' ]]; then
        brew install $1
    fi

    if is_installed $1; then
        echo "$1 installed successfully"
    else
        echo "Was not able to install $1 automatically. Please install manually"
    fi
}

unamestr=$(uname)

# Try to install the external programs automatically
install_program mmseqs
install_program mafft
install_program fasttree
install_program iqtree