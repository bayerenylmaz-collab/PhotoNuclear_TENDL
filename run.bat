@echo off
setlocal EnableExtensions
cd /d "%~dp0"

echo.
echo ============================================
echo  Photonuclear Katalog - kolay calistirici
echo ============================================
echo.

rem Find a real Python (not the Microsoft Store stub).
set "PYEXE="
set "PYARGS="

where py >nul 2>nul
if not errorlevel 1 (
  py -3 -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
  if not errorlevel 1 (
    set "PYEXE=py"
    set "PYARGS=-3"
    goto :python_ok
  )
)

where python >nul 2>nul
if not errorlevel 1 (
  python -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
  if not errorlevel 1 (
    set "PYEXE=python"
    set "PYARGS="
    goto :python_ok
  )
)

where python3 >nul 2>nul
if not errorlevel 1 (
  python3 -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
  if not errorlevel 1 (
    set "PYEXE=python3"
    set "PYARGS="
    goto :python_ok
  )
)

echo [HATA] Calisan bir Python 3.10+ bulunamadi.
echo.
echo Windows'ta sik gorulen neden: Python yok veya Microsoft Store
echo kisayolu acik, gercek Python kurulu degil.
echo.
echo Yapman gerekenler:
echo   1) https://www.python.org/downloads/ adresinden Python 3.12/3.13 indir
echo   2) Kurulumda mutlaka "Add python.exe to PATH" isaretle
echo   3) Bu pencereyi kapat, bilgisayari yeniden baslatmana gerek yok;
echo      yeni bir run.bat penceresi acman yeterli
echo.
echo Alternatif: WSL kullanip ./run.sh calistirabilirsin.
echo.
pause
exit /b 1

:python_ok
echo Python bulundu: %PYEXE% %PYARGS%
%PYEXE% %PYARGS% -c "import sys; print('  surum:', sys.version.split()[0])"
echo.

if not exist ".venv\Scripts\python.exe" (
  echo Sanal ortam olusturuluyor...
  %PYEXE% %PYARGS% -m venv .venv
  if errorlevel 1 (
    echo [HATA] venv olusturulamadi.
    echo Python kurulumunu kontrol et ^(PATH + pip/venv^).
    pause
    exit /b 1
  )
)

call ".venv\Scripts\activate.bat"
python -m pip install --upgrade pip >nul
echo Paket kuruluyor / guncelleniyor...
python -m pip install -e .
if errorlevel 1 (
  echo [HATA] paket kurulumu basarisiz.
  pause
  exit /b 1
)

echo.
python -m photonuclear_catalog --interactive
set ERR=%ERRORLEVEL%
echo.
if not "%ERR%"=="0" (
  echo [HATA] program %ERR% kodu ile bitti.
) else (
  echo Tamam. Cikti klasorundeki report.html dosyasini acabilirsiniz.
)
echo.
pause
exit /b %ERR%
