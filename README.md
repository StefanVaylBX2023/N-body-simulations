# N-body-simulations

To perform visualization of the results the SFML library is used.
To run the project make sure you have installed SFML library. If not, please, proceed with installation with homebrew (macOS) or using this link \href{https://www.sfml-dev.org/download.php}{link}. \\
For the compilation you will need to specify dependencies of the library to the compiler:
\begin{enumerate}
    \item Find the installation path of SFML (for example if you used brew: "brew info sfml"
    \item explicitly find locations of "lib" and "include" directory
\end{enumerate}
To perform compilation one may need to specify the following flags for the compiler:
\begin{enumerate}
    \item  \texttt{ main.cpp -o main}
    \item \texttt{-pthread} (to enable multithreading)
    \item \texttt{-std=c++11}
    \item \texttt{-I/path/to/sfml/include} (insert your path to "include")
    \item \texttt{-L/path/to/sfml/lib}  (insert your path to "lib")
    \item \texttt{-lsfml-graphics -lsfml-window -lsfml-system} (linking of library dependencies)
\end{enumerate}
The final compilation prompt then look the following way:\\
\texttt{g++ main.cpp -o main -std=c++11 -pthread -I/path/to/sfml/include\\ -L/path/to/sfml/lib -lsfml-graphics -lsfml-window -lsfml-system}\\
Example: \\
\texttt{g++ main.cpp -I/opt/homebrew/Cellar/sfml/2.5.1\_2/include\\ -L/opt/homebrew/Cellar/sfml/2.5.1\_2/lib -o main -lsfml-graphics \\-lsfml-window -lsfml-system -std=c++11 -pthread}
