$RepoRoot  = $PSScriptRoot
$ResultDirs = Get-ChildItem -Path $RepoRoot -Directory -Filter "Result*"

if ($ResultDirs.Count -eq 0) {
    Write-Error "No Result* directories found in $RepoRoot"
    exit 1
}

Write-Host "Found $($ResultDirs.Count) Result* directories to process"
Write-Host "Skipping any subfolder starting with 'Analysis'"
Write-Host ""

$totalCompressed = 0
$totalErrors     = 0

foreach ($resultDir in $ResultDirs) {
    Write-Host "==> $($resultDir.Name)"

    $folders = Get-ChildItem -Path $resultDir.FullName -Directory |
               Where-Object { $_.Name -notlike "Analysis*" }

    foreach ($folder in $folders) {
        $output = Join-Path $resultDir.FullName "$($folder.Name).tar.gz"
        Write-Host "  Compressing: $($folder.Name)"

        tar -czf $output -C $resultDir.FullName $folder.Name

        if ($LASTEXITCODE -eq 0) {
            $size = [math]::Round((Get-Item $output).Length / 1MB, 2)
            Write-Host "    Done ($size MB)"
            $totalCompressed++
        } else {
            Write-Host "    ERROR compressing $($folder.Name)"
            $totalErrors++
        }
    }
    Write-Host ""
}

Write-Host "Finished. Compressed: $totalCompressed  Errors: $totalErrors"
