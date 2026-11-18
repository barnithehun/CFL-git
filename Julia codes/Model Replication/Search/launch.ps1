# launch.ps1 <wave> [<workers>]: starts one Julia process per worker, each running its share of the wave's candidates
# (candidates.csv rows whose wave column equals <wave>, dealt out in order), through Search.jl run <ids>.
# Output of worker i goes to logs/<wave>_worker<i>.txt; each candidate also writes runs/<id>/log.txt and DIAG.csv.
param([Parameter(Mandatory = $true)][string]$wave, [int]$workers = 8)

$search = $PSScriptRoot
$repo   = (Resolve-Path (Join-Path $search "..\..\..")).Path
$julia  = "C:\Users\szjud\.julia\juliaup\julia-1.12.7+0.x64.w64.mingw32\bin\julia.exe"
$logs   = Join-Path $search "logs"
New-Item -ItemType Directory -Force $logs | Out-Null

$rows = @(Import-Csv (Join-Path $search "candidates.csv") | Where-Object { $_.wave -eq $wave })
if ($rows.Count -eq 0) { Write-Host "no candidates in wave $wave"; exit 1 }
$ids = @($rows | ForEach-Object { $_.id })
Write-Host "Wave $wave : $($ids.Count) candidates on $workers workers: $($ids -join ', ')"

for ($i = 0; $i -lt $workers; $i++) {
    $mine = @()
    for ($j = $i; $j -lt $ids.Count; $j += $workers) { $mine += $ids[$j] }
    if ($mine.Count -eq 0) { continue }
    $log = Join-Path $logs "${wave}_worker$i.txt"
    $err = Join-Path $logs "${wave}_worker$i.err.txt"
    $args = @("--project=`"$repo`"", "`"$search\Search.jl`"", "run") + $mine
    Write-Host "worker $i : $($mine -join ', ')  -> $log"
    Start-Process -FilePath $julia -ArgumentList $args -NoNewWindow -RedirectStandardOutput $log -RedirectStandardError $err
    Start-Sleep -Seconds 3
}
Write-Host "started; the Julia processes keep running after this script ends"
