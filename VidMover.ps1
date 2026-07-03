$SourceDir = "E:\GE\GE Ethanol Quadrants"
$DestDir = "E:\GE\GE Ethanol Quadrants resampled"
$TargetSuffix = "*output_sampled.avi"

if (-Not (Test-Path $SourceDir)) {
    Write-Host "Error: Source directory '$SourceDir' does not exist." -ForegroundColor Red
    exit
}

if (-Not (Test-Path $DestDir)) {
    New-Item -ItemType Directory -Path $DestDir -Force | Out-Null
}

Get-ChildItem -Path $SourceDir -Directory -Recurse | ForEach-Object {
    $relativePath = $_.FullName.Substring($SourceDir.Length + 1)
    $targetFolder = Join-Path -Path $DestDir -ChildPath $relativePath
    if (-Not (Test-Path $targetFolder)) {
        New-Item -ItemType Directory -Path $targetFolder -Force | Out-Null
    }
}

Write-Host "Moving '$TargetSuffix' files..." -ForegroundColor Cyan

$filesToMove = Get-ChildItem -Path $SourceDir -Filter $TargetSuffix -File -Recurse
$movedCount = 0

foreach ($file in $filesToMove) {
    $relativePath = $file.FullName.Substring($SourceDir.Length + 1)
    $targetFilePath = Join-Path -Path $DestDir -ChildPath $relativePath
    Move-Item -Path $file.FullName -Destination $targetFilePath -Force
    $movedCount++
}

Write-Host "all done!" -ForegroundColor Green
