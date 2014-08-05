Pairwise Asymmetric Interference 
===

This code performs the convergent-cross mapping (CCM) and pairwise asymmetric interference (PAI) calculations in "Convergent Cross-Mapping and Pairwise Asymmetric Inference", James M. McCracken and Robert S. Weigel [http://arxiv.org/abs/1407.5696] (doi: [TBD]).

The code calculates CCM correlations slightly differently than Sugihara et al., who originally introduced the techniques in "Detecting Causality in Complex Ecosystems" [http://www.sciencemag.org/content/338/6106/496.figures-only].  This code uses the square of the correlation between the estimated time series and the original time series, rather than the unsquared version used by Sugihara et al.  See the above paper for details on the algorithms used in this code.    

This code requires C++11. To compile, enter

    make all

To check the compile, enter

    make test

To run, enter

    ./PAI -h

to see the command line options and requirements.  All of the input and output files are currently plain text with the specific formats indicated by PAI -h.

This code has been compiled and tested on Ubuntu 13.10 with kernel version 3.11.0-26-generic and GCC with the following version information:
   
    $ gcc -v
    Using built-in specs.
    COLLECT_GCC=gcc
    COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-linux-gnu/4.8/lto-wrapper
    Target: x86_64-linux-gnu
    Configured with: ../src/configure -v --with-pkgversion='Ubuntu/Linaro 4.8.1-10ubuntu9' --with-bugurl=file:///usr/share/doc/gcc-4.8/README.Bugs --enable-languages=c,c++,java,go,d,fortran,objc,obj-c++ --prefix=/usr --program-suffix=-4.8 --enable-shared --enable-linker-build-id --libexecdir=/usr/lib --without-included-gettext --enable-threads=posix --with-gxx-include-dir=/usr/include/c++/4.8 --libdir=/usr/lib --enable-nls --with-sysroot=/ --enable-clocale=gnu --enable-libstdcxx-debug --enable-libstdcxx-time=yes --enable-gnu-unique-object --enable-plugin --with-system-zlib --disable-browser-plugin --enable-java-awt=gtk --enable-gtk-cairo --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-4.8-amd64/jre --enable-java-home --with-jvm-root-dir=/usr/lib/jvm/java-1.5.0-gcj-4.8-amd64 --with-jvm-jar-dir=/usr/lib/jvm-exports/java-1.5.0-gcj-4.8-amd64 --with-arch-directory=amd64 --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --enable-objc-gc --enable-multiarch --disable-werror --with-arch-32=i686 --with-abi=m64 --with-multilib-list=m32,m64,mx32 --with-tune=generic --enable-checking=release --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu
    Thread model: posix
    gcc version 4.8.1 (Ubuntu/Linaro 4.8.1-10ubuntu9)
