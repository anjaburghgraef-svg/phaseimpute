@echo off
REM Windows batch file to launch the Shiny app

echo ========================================
echo phaseimpute QC Report
echo ========================================
echo.

REM Check if data directory is provided as argument
set DATA_DIR=%1
if "%DATA_DIR%"=="" set DATA_DIR=run_testprofile/beagle5

REM Check if port is provided as argument
set PORT=%2
if "%PORT%"=="" set PORT=3838

echo Data directory: %DATA_DIR%
echo Port: %PORT%
echo URL: http://127.0.0.1:%PORT%
echo.
echo Press Ctrl+C to stop the server
echo ========================================
echo.

REM Launch the app using full path to Rscript
"C:\Program Files\R\R-4.4.2\bin\Rscript.exe" launch_report.R "%DATA_DIR%" %PORT%

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ========================================
    echo Error launching app
    echo Please check that dependencies are installed:
    echo   install_dependencies.bat
    echo ========================================
    pause
)
