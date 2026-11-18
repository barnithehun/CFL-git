# make_target_variants.ps1: runs calib_targets.R once per target variant used by the calibration search
# (Julia codes/Model Replication/Search/candidates.csv, column targets). Each run writes calib_targets_<tag>.csv and
# calib_targets_<tag>_settings.csv next to calib_targets.csv; the default file is not touched.
$rscript = "C:\Program Files\R\R-4.5.2\bin\Rscript.exe"
$here = $PSScriptRoot
$variants = @(
    @("A8",    "size_var=assets", "cutoff_usd=1e8"),
    @("L9",    "size_var=ltq",    "cutoff_usd=1e9"),
    @("L7",    "size_var=ltq",    "cutoff_usd=1e7"),
    @("L8m6",  "size_var=ltq",    "cutoff_usd=1e8", "micro_usd=1e6"),
    @("L8d2",  "size_var=ltq",    "cutoff_usd=1e8", "d2c_cut=2"),
    @("L8d4",  "size_var=ltq",    "cutoff_usd=1e8", "d2c_cut=4"),
    @("L8y20", "size_var=ltq",    "cutoff_usd=1e8", "end_year=2020"),
    @("A9",    "size_var=assets", "cutoff_usd=1e9")
)
foreach ($v in $variants) {
    $tag = $v[0]
    $args = @("calib_targets.R") + $v[1..($v.Count - 1)] + @("out_tag=$tag")
    Write-Host "=== $tag : $($args -join ' ')"
    & $rscript $args 2>&1 | Select-Object -Last 12
}
Write-Host "done"
