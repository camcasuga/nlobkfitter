#nlobkfitter installation

## Dependencies:
    - cmake
    - GSL
    - amplitudelib_v2
    - Minuit2
    - Cuba (MC integration library)

## Procedure

    1. install GSL
    2. build amplitudelib_v2 
       - clone https://github.com/hejajama/amplitudelib
       - mkdir [amplitudelib/]build && cd build && cmake .. && make
    3. build Minuit2:
       - download from CERN homepage
       - ./configure --disable-openmp --prefix=INSTALLDIR
       - make -j[N] && make install
    4. build Cuba:
       - download
       - ./configure --prefix=INSTALLDIR
       - make && make install
    5. build nlobkfitter & tools:
       - clone https://github.com/hejajama/nlobkfitter
       - cd build && make
