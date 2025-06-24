# Seismic Wave Propagation Simulator

## Descripción del Proyecto

Este proyecto simula la propagación de ondas sísmicas en un medio 2D utilizando la ecuación de onda. Una fuente de energía, modelada como una wavelet de Ricker, excita el medio en un punto específico. La simulación calcula la amplitud de la onda en cada punto de una malla a lo largo del tiempo.

El programa principal está escrito en C++ y utiliza CMake para la compilación. Los parámetros de la simulación se configuran a través de un archivo `config.json`. La salida de la simulación consiste en una serie de archivos de datos que representan el estado de la malla en diferentes instantes de tiempo.

Se incluye un script de Python (`visualize.py`) para procesar estos archivos de datos y generar una animación en formato GIF que muestra la propagación de las ondas en 3D.

## Dependencias

### Para Compilar (C++)
*   **CMake** (versión 3.10 o superior)
*   Un **compilador de C++** compatible con C++17 (ej. GCC, Clang, MSVC)
*   (Opcional, pero recomendado para `main.cpp`) Una **librería para parsear JSON** como [nlohmann/json](https://github.com/nlohmann/json). El esqueleto actual de `main.cpp` simula la lectura de JSON, pero para una funcionalidad completa, se debería integrar una librería de este tipo. Para usar nlohmann/json, puedes incluir su archivo `json.hpp` directamente en tu proyecto o enlazarla mediante CMake.

### Para Visualizar (Python)
*   **Python** (versión 3.6 o superior)
*   Las siguientes librerías de Python (instalar con `pip install -r requirements.txt`):
    *   `numpy`
    *   `pandas`
    *   `matplotlib`
    *   `imageio` (y potencialmente `imageio[ffmpeg]` para una mejor compatibilidad con GIFs)

## Cómo Compilar (C++)

1.  **Asegúrate de tener CMake y un compilador de C++ instalados.**
2.  **Si usas una librería externa para JSON (como nlohmann/json):**
    *   Descarga `json.hpp` y colócala en un directorio accesible (ej. `src/include/nlohmann/json.hpp`).
    *   O sigue las instrucciones de la librería para integrarla con CMake (ej. usando `find_package` si está instalada en el sistema, o `add_subdirectory` si la incluyes como un submodulo git).
    *   Deberás ajustar el `CMakeLists.txt` para que encuentre y enlace la librería JSON. El `main.cpp` actual tiene comentarios sobre esto.
3.  **Crea un directorio de compilación y navega hacia él:**
    ```bash
    mkdir build
    cd build
    ```
4.  **Ejecuta CMake para configurar el proyecto:**
    Asegúrate de que el `CMakeLists.txt` raíz esté configurado para encontrar tu código fuente (ej. `src/main.cpp`) y cualquier librería necesaria.
    ```bash
    cmake ..
    ```
    *Nota: Si has colocado `nlohmann/json.hpp` en `src/include`, y tu `main.cpp` la incluye como `#include "nlohmann/json.hpp"`, puede que necesites añadir `target_include_directories(seismic_sim PUBLIC ../src/include)` al `CMakeLists.txt` (después de `add_executable`). El esqueleto del `CMakeLists.txt` tiene comentarios al respecto.*

5.  **Compila el proyecto:**
    ```bash
    make
    ```
    (O el comando equivalente para tu sistema de compilación, ej. `ninja` o abrir el proyecto en un IDE como Visual Studio).

    Después de una compilación exitosa, deberías encontrar el ejecutable (ej. `seismic_sim`) en el directorio `build` o `build/src`.

## Cómo Ejecutar

### 1. Ejecutar la Simulación C++

*   Una vez compilado, el ejecutable (por ejemplo, `seismic_sim`) se encontrará típicamente en el directorio `build`.
*   Asegúrate de que el archivo `config.json` esté presente en el directorio raíz del proyecto (el mismo nivel que `src`, `data`, `build`).
*   Desde el directorio raíz del proyecto, puedes ejecutar la simulación:
    ```bash
    ./build/seismic_sim
    ```
    (Ajusta la ruta al ejecutable si es necesario, por ejemplo, en Windows podría ser `build\Debug\seismic_sim.exe`).
*   La simulación leerá `config.json` y guardará los archivos de datos (ej. `frame_0000.dat`, `frame_0001.dat`, etc.) en el directorio `/data`.

### 2. Generar la Animación (Python)

*   Asegúrate de haber instalado las dependencias de Python listadas en `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```
*   Una vez que la simulación C++ haya generado los archivos de datos en el directorio `/data`, ejecuta el script de visualización desde el directorio raíz del proyecto:
    ```bash
    python visualize.py
    ```
*   El script buscará archivos en `/data`, procesará cada uno para generar un cuadro de la animación, y finalmente guardará la animación como `simulation_output.gif` en el directorio raíz.
*   Si el directorio `/data` está vacío o no contiene los archivos esperados, el script mostrará un mensaje de error.

---
*Este es un esqueleto inicial. Necesitarás implementar la lógica completa de la simulación en `src/main.cpp` y potencialmente refinar el script `visualize.py` y el `CMakeLists.txt` para una funcionalidad completa y robusta.*
