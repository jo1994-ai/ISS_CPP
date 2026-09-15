@echo off
rem Build and run the Q15 FFT example on Windows.
rem Needs Visual Studio 2019/2022 with "Desktop development with C++".
setlocal
cd /d "%~dp0"

set "CMAKE=cmake"
where cmake >nul 2>nul
if %errorlevel%==0 goto have_cmake
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" goto no_vs
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -find Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe`) do set "CMAKE=%%i"
if "%CMAKE%"=="cmake" goto no_vs

:have_cmake
if not exist build\CMakeCache.txt (
  "%CMAKE%" -S . -B build
  if errorlevel 1 goto failed
)
"%CMAKE%" --build build --config Release
if errorlevel 1 goto failed

echo.
fft_example_q15.exe %*
if errorlevel 1 goto failed
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
