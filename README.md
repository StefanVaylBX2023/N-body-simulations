# N-body-simulations

o perform visualization of the results, the SFML library is used. To run the project, make sure you have installed the SFML library. If not, please proceed with the installation using Homebrew (macOS) or by following this link: SFML Download.

For compilation, you will need to specify the dependencies of the library to the compiler. Follow these steps:

Find the installation path of SFML. For example, if you used Homebrew, run the command: brew info sfml.
Explicitly find the locations of the "lib" and "include" directories.
To compile the project, you may need to specify the following flags for the compiler:

main.cpp -o main
-pthread (to enable multithreading)
-std=c++11
-I/path/to/sfml/include (insert your path to the "include" directory)
-L/path/to/sfml/lib (insert your path to the "lib" directory)
-lsfml-graphics -lsfml-window -lsfml-system (linking of library dependencies)

