# Toolchain for the bundled w64devkit (MinGW-w64 GCC, msvcrt) on Windows.
#
# It also bypasses CMake's compiler-probe step.  On this machine Windows
# Defender briefly locks freshly produced executables, which makes CMake's
# "read a.exe to detect the compiler" step fail intermittently
# ("file failed to open for reading (Invalid argument)").  Presetting the
# compiler identity and marking the compilers as FORCED skips that probe.

set(CMAKE_SYSTEM_NAME Windows)

set(W64DEVKIT "D:/Program Files/MinGW/w64devkit/bin")
set(CMAKE_C_COMPILER   "${W64DEVKIT}/gcc.exe")
set(CMAKE_CXX_COMPILER "${W64DEVKIT}/g++.exe")
set(CMAKE_MAKE_PROGRAM "${W64DEVKIT}/mingw32-make.exe")

# --- skip compiler identification / ABI probing (Defender race) ---
set(CMAKE_C_COMPILER_ID GNU)
set(CMAKE_CXX_COMPILER_ID GNU)
set(CMAKE_C_COMPILER_VERSION 14.2.0)
set(CMAKE_CXX_COMPILER_VERSION 14.2.0)
set(CMAKE_C_COMPILER_ID_RUN TRUE)
set(CMAKE_CXX_COMPILER_ID_RUN TRUE)
set(CMAKE_C_COMPILER_FORCED TRUE)
set(CMAKE_CXX_COMPILER_FORCED TRUE)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_SIZEOF_VOID_P 8)
set(CMAKE_C_SIZEOF_DATA_PTR 8)
set(CMAKE_CXX_SIZEOF_DATA_PTR 8)
