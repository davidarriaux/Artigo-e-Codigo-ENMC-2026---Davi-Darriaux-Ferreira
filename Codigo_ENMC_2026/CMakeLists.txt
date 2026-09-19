cmake_minimum_required(VERSION 3.10)
project(SOR_PNBI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# O Eigen (biblioteca de algebra linear) ja vem dentro de src/Eigen: nao e preciso instalar nada.
add_executable(sor_pnbi
        src/main.cpp
        src/Problema.cpp
        src/NBI.cpp
        src/Relaxacao.cpp
        src/Metodos.cpp
)

# Permite ao programa achar a pasta Problemas/ mesmo rodando a partir da pasta de build
target_compile_definitions(sor_pnbi PRIVATE RAIZ_PROJETO="${CMAKE_SOURCE_DIR}")

# Otimizacao maxima (tambem no perfil Debug do CLion, senao o programa fica muito lento)
if(MSVC)
    target_compile_options(sor_pnbi PRIVATE /O2 /utf-8)
else()
    target_compile_options(sor_pnbi PRIVATE -O3)
endif()
