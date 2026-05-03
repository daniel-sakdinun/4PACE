@echo off
echo 🚀 Starting 4PACE Release Process...

:: 1. Clean up old build artifacts
echo 🧹 Cleaning old builds...
if exist dist ( rmdir /s /q dist )
if exist build ( rmdir /s /q build )
del /s /q *.egg-info >nul 2>&1

:: 2. Build the package
echo 📦 Building package (sdist and wheel)...
python -m build

:: 3. Upload to PyPI
echo 📤 Uploading to PyPI via Twine...
python -m twine upload dist/*

echo.
echo ✅ Done! Your new version of 4PACE is now live on PyPI.
pause