@echo off
rem Install the PRESTUS Nexus (Windows) in one step: double-click this file,
rem or run  nexus\install.bat  in a command prompt.
rem   nexus\install.bat /y    install Miniforge without asking and do not pause
rem
rem Uses conda/mamba if it is installed; otherwise offers to install Miniforge
rem (conda-forge's conda) into %USERPROFILE%\miniforge3. Then creates the
rem environment 'prestus-nexus' from environment.yml, or updates it if it
rem exists. Run it again after a PRESTUS update to update the environment.
setlocal EnableExtensions
set "HERE=%~dp0"
set "ENV_NAME=prestus-nexus"
set "CONDA="
set "YES="
if /i "%~1"=="/y" set "YES=1"
if /i "%~1"=="--yes" set "YES=1"

for %%C in (mamba.bat conda.bat mamba.exe conda.exe) do (
    if not defined CONDA for /f "delims=" %%P in ('where %%C 2^>nul') do if not defined CONDA set "CONDA=%%P"
)
for %%D in (miniforge3 mambaforge miniconda3 anaconda3) do (
    if not defined CONDA if exist "%USERPROFILE%\%%D\condabin\conda.bat" set "CONDA=%USERPROFILE%\%%D\condabin\conda.bat"
)
if not defined CONDA if exist "%LOCALAPPDATA%\miniforge3\condabin\conda.bat" set "CONDA=%LOCALAPPDATA%\miniforge3\condabin\conda.bat"
if not defined CONDA if exist "%LOCALAPPDATA%\anaconda3\condabin\conda.bat" set "CONDA=%LOCALAPPDATA%\anaconda3\condabin\conda.bat"
if defined CONDA goto :have_conda

echo conda was not found. The Nexus needs it once, to install its Python environment.
if defined YES goto :install_miniforge
choice /m "Install Miniforge into %USERPROFILE%\miniforge3 now"
if errorlevel 2 goto :declined

:install_miniforge
set "MF_EXE=%TEMP%\Miniforge3-Windows-x86_64.exe"
echo Downloading Miniforge...
curl -fsSL -o "%MF_EXE%" https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe
if errorlevel 1 goto :fail
echo Installing Miniforge: this takes a few minutes and opens no window.
start /wait "" "%MF_EXE%" /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%USERPROFILE%\miniforge3
del "%MF_EXE%"
set "CONDA=%USERPROFILE%\miniforge3\condabin\conda.bat"
if not exist "%CONDA%" (
    echo Miniforge did not install to %USERPROFILE%\miniforge3.
    goto :fail
)

:have_conda
echo Using %CONDA%
pushd "%HERE%"
call "%CONDA%" env list | findstr /b /c:"%ENV_NAME% " >nul
if errorlevel 1 goto :create
echo Updating the environment %ENV_NAME%...
call "%CONDA%" env update -n %ENV_NAME% -f environment.yml
goto :created
:create
echo Creating the environment %ENV_NAME%...
call "%CONDA%" env create -n %ENV_NAME% -f environment.yml
:created
if errorlevel 1 goto :fail_popd

call "%CONDA%" run -n %ENV_NAME% python -c "import yaml, PyQt6.QtWidgets, prestus_nexus"
if errorlevel 1 goto :fail_import
for /f "delims=" %%P in ('call "%CONDA%" run -n %ENV_NAME% python -c "import sys; print(sys.prefix)"') do set "ENV_DIR=%%P"
popd
echo.
echo PRESTUS Nexus installed.
echo   Start it in MATLAB:  prestus_nexus
echo   or directly:         "%ENV_DIR%\Scripts\prestus-nexus.exe"
if not defined YES pause
exit /b 0

:declined
echo Stopped. Install conda yourself, then run this again.
goto :fail

:fail_import
echo The environment was created, but it cannot import the Nexus.
echo See Installation in doc\doc_nexus.md.
:fail_popd
popd
:fail
if not defined YES pause
exit /b 1
