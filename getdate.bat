@echo off
setlocal EnableDelayedExpansion

:: Extract the date into YYYY MM DD variables
for /f "tokens=2-4 delims=/ " %%a in ('echo %date%') do (
    set year=%%c
    set month=%%a
    set day=%%b
)

:: Add leading zeros to month and day
if %month% lss 10 set month=0%month%
if %day% lss 10 set day=0%day%

:: Combine to form YYYYMMDD
set formattedDate=%year%%month%%day%

echo date=%formattedDate%

endlocal
