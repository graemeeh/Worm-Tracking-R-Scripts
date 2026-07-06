# moves all items of a given extension (for me, all .avi files) to recycle bin recursively
# in current folder and all subfolders


Add-Type -AssemblyName Microsoft.VisualBasic
if (-NOT ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole] "Administrator"))  
{  
  $arguments = "& '" +$myinvocation.mycommand.definition + "'"
  Start-Process powershell -Verb runAs -ArgumentList $arguments
  Break
}


function Remove-Item-toRecycle($item) {
    Get-Item -Path $item | %{ $fullpath = $_.FullName}
    [Microsoft.VisualBasic.FileIO.FileSystem]::DeleteFile($fullpath,'OnlyErrorDialogs','SendToRecycleBin')
}


Get-ChildItem -Path $PSScriptRoot *.avi -Recurse | foreach { Remove-Item-toRecycle -item $_.FullName }
foreach ($s in Get-ChildItem -Path $PSScriptRoot *.avi -Recurse)
{
    Write-Host "s: " $s
}