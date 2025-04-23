@echo off
echo ===========================
echo Building the project...
echo ===========================
cmake --build build

echo.
echo ===========================
echo Running heartrate.exe...
echo ===========================
build\Debug\heartrate.exe INPUT_Dataset OUTPUT_Dataset

echo.
echo ===========================
echo Execution finished.
echo Press any key to close.
echo ===========================
pause
