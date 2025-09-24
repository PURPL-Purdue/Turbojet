echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v252\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\tell.exe" DESKTOP-LG4V0VL 64320 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 26664) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 26644) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 7808) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 20988) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 22556) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 26180)
del "C:\Users\rjmcm\Documents\GitHub\Turbojet\Combustor\Simulation\9_23_25_Axisymmetric_Full_Combustor\cleanup-fluent-DESKTOP-LG4V0VL-22556.bat"
