#!/bin/bash

function format {
    echo "Formatting cpp files in ${1}"
    find "$1" -name "*.cpp" -exec clang-format -i -style=file {} \;
    find "$1" -name "*.h" -exec clang-format -i -style=file {} \;
}

format "."
