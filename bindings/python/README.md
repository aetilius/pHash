# Python binding to use pHash on images

## Local use

Check that all dependencies are available and create the header pHash.h:

```
mkdir build-phash
cd build-phash
cmake ../../../
cd ..
```

Build the python extension in the current directory:

```
python setup.py build_ext --inplace
```


