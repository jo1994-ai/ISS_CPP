@echo off
rem Build and run fft_example on Windows.
rem Needs Visual Studio 2019 or 2022 with "Desktop development with C++"
rem (that includes the C++ compiler and CMake). FFTW is bundled.
setlocal
cd /d "%~dp0"

rem ---- find cmake: on PATH, or the copy that ships with Visual Studio
set "CMAKE=cmake"
where cmake >nul 2>nul
if %errorlevel%==0 goto have_cmake

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" goto no_vs
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -find Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe`) do set "CMAKE=%%i"
if "%CMAKE%"=="cmake" goto no_vs

:have_cmake
echo Using CMake: %CMAKE%

rem ---- configure (first time: extracts and prepares FFTW)
if not exist build\CMakeCache.txt (
  "%CMAKE%" -S . -B build
  if errorlevel 1 goto failed
)

rem ---- build (first time: compiles FFTW, a few minutes)
"%CMAKE%" --build build --config Release --parallel
if errorlevel 1 goto failed

rem ---- run
echo.
fft_example.exe %*
if errorlevel 1 goto failed
echo.
echo Done. Reconstructed audio: output\reconstructed.wav
goto end

:no_vs
echo ERROR: CMake / Visual Studio not found.
echo Install Visual Studio 2022 Community with "Desktop development with C++":
echo   https://visualstudio.microsoft.com/downloads/
goto end_error

:failed
echo.
echo ERROR: build or run failed, see the messages above.

:end_error
pause
exit /b 1

:end
pause
