<#
  Minimal sweep: ONLY set EHEP_BOTTOM_LEN_M
  Run:
    powershell -NoProfile -ExecutionPolicy Bypass -File .\sweep.ps1
#>
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$ScriptDir = if ($PSScriptRoot) { $PSScriptRoot } else { Split-Path -Parent $MyInvocation.MyCommand.Path }

function Find-EnhanceExe([string]$baseDir){
  if ($env:ENHANCE_EXE -and (Test-Path $env:ENHANCE_EXE)) { return (Resolve-Path $env:ENHANCE_EXE).Path }
  $candidates = @(
    (Join-Path $baseDir 'x64\Release\enhance.exe'),
    (Join-Path $baseDir 'Release\enhance.exe'),
    (Join-Path $baseDir 'enhance.exe'),
    (Join-Path $baseDir 'x64\Debug\enhance.exe'),
    (Join-Path $baseDir 'Debug\enhance.exe')
  )
  foreach($p in $candidates){ if(Test-Path $p){ return (Resolve-Path $p).Path } }
  $found = Get-ChildItem -Path $baseDir -Filter 'enhance.exe' -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
  if($found){ return $found.FullName }
  throw 'enhance.exe not found. Build Release|x64 or set ENHANCE_EXE.'
}

$exe = Find-EnhanceExe $ScriptDir
Write-Host "[Info] exe: $exe"

# ONLY sweep this
$enhBottomLenList = @(200,400,600,800,1000,1200)

$ts = Get-Date -Format 'yyyyMMdd_HHmmss'
$runRoot = Join-Path $ScriptDir ("runs_"+$ts)
New-Item -ItemType Directory -Force -Path $runRoot | Out-Null

# one-line-per-case summary
$summaryFile = Join-Path $runRoot 'summary_runs.csv'
"enhB_m,Q_load_kWh,Q_src_kWh,P_el_kWh,P_pump_kWh,COP_annual,HP_on_hours,results_csv" | Set-Content -Path $summaryFile -Encoding UTF8

function ToDouble($v){
  try { $d = [double]$v } catch { return 0.0 }
  if([double]::IsNaN($d) -or [double]::IsInfinity($d)){ return 0.0 }
  return $d
}

function Run-OneCase([double]$enhB){
  # THE ONLY ENV VAR YOU SET:
  $env:EHEP_BOTTOM_LEN_M = [string]$enhB

  Write-Host "[Run] EHEP_BOTTOM_LEN_M=$enhB"
  Push-Location $ScriptDir
    & $exe
  Pop-Location
  if($LASTEXITCODE -ne 0){ Write-Warning "process exit code $LASTEXITCODE" }

  # save results
  $tag = "enhB_${enhB}m"
  $resOut = Join-Path $runRoot ("results_"+$tag+".csv")
  if(Test-Path (Join-Path $ScriptDir 'results.csv')){
    Copy-Item (Join-Path $ScriptDir 'results.csv') $resOut -Force
  }

  # summarize (simple sum of kW columns over rows; assumes results is hourly or already kWh-compatible per row)
  if(Test-Path $resOut){
    $rows = Import-Csv -Path $resOut
    if($rows.Count -gt 0){
      $sumQspace=0.0; $sumQdhw=0.0; $sumPel=0.0; $sumPpump=0.0; $sumQsrc=0.0; $onRows=0

      foreach($r in $rows){
        if($r.PSObject.Properties.Name -contains 'Q_space_served_kW'){ $sumQspace += ToDouble $r.Q_space_served_kW }
        if($r.PSObject.Properties.Name -contains 'Q_dhw_served_kW'){   $sumQdhw   += ToDouble $r.Q_dhw_served_kW }
        if($r.PSObject.Properties.Name -contains 'P_el_kW'){            $sumPel    += ToDouble $r.P_el_kW }
        if($r.PSObject.Properties.Name -contains 'P_pump_kW'){          $sumPpump  += ToDouble $r.P_pump_kW }
        if($r.PSObject.Properties.Name -contains 'Q_geo_kW'){           $sumQsrc   += ToDouble $r.Q_geo_kW }
        if($r.PSObject.Properties.Name -contains 'HP_on'){
          $on = 0; try { $on = [int](ToDouble $r.HP_on) } catch {}
          if($on -gt 0){ $onRows += 1 }
        }
      }

      $Q_load_kWh = $sumQspace + $sumQdhw
      $Q_src_kWh  = $sumQsrc
      $P_el_kWh   = $sumPel
      $P_pump_kWh = $sumPpump
      $COP_annual = if(($P_el_kWh + $P_pump_kWh) -gt 1e-9){ [double]($Q_load_kWh / ($P_el_kWh + $P_pump_kWh)) } else { 0.0 }
      $HP_on_hours = [double]$onRows

      Add-Content -Path $summaryFile -Value ("$enhB,$Q_load_kWh,$Q_src_kWh,$P_el_kWh,$P_pump_kWh,$COP_annual,$HP_on_hours,$resOut")
    }
  }
}

foreach($enhB in $enhBottomLenList){
  Run-OneCase -enhB $enhB
}

Write-Host "[Done] output: $runRoot"
