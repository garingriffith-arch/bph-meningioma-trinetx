# Run from this repository root on Windows. Requires git, gh, and Python.
$ErrorActionPreference = 'Stop'
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) { throw 'Install GitHub CLI (gh) first.' }
if (-not (Get-Command git -ErrorAction SilentlyContinue)) { throw 'Install Git first.' }
& gh auth status
if ($LASTEXITCODE -ne 0) { throw 'Run gh auth login.' }
if (-not (Test-Path 'R/27_prepare_resident_release.R')) { throw 'Run from the repository root.' }
if (Get-Command python -ErrorAction SilentlyContinue) { & python tests/check_repository.py; if ($LASTEXITCODE -ne 0) { throw 'Repository hygiene check failed.' } }
$repo = 'garingriffith-arch/bph-meningioma-trinetx'
& gh repo view $repo --json name,url,isPrivate 2>$null
$exists = $LASTEXITCODE -eq 0
if ($exists) { throw "Remote $repo already exists. Review its history before merging or pushing; no automatic overwrite." }
if (-not (Test-Path .git)) { & git init -b main; if ($LASTEXITCODE -ne 0) { throw 'git init failed.' } }
& git add -- README.md .gitattributes .gitignore .github R docs reporting source-provenance scripts tests
if ($LASTEXITCODE -ne 0) { throw 'git add failed.' }
Write-Host 'Files to publish to a PRIVATE repository:'
& git diff --cached --name-only
Write-Host 'This repository must NOT include restricted data, output ZIPs, or patient-level files.'
$confirm = Read-Host 'Enter PRIVATE to confirm the staged source files have been reviewed'
if ($confirm -cne 'PRIVATE') { throw 'Stopped without pushing.' }
& git commit -m 'Document reproducible meningioma pharmacoepidemiology analysis'
if ($LASTEXITCODE -ne 0) { throw 'Commit failed.' }
& gh repo create $repo --private --source . --remote origin --push
if ($LASTEXITCODE -ne 0) { throw 'GitHub repository creation or push failed.' }
& gh repo view $repo --json url,isPrivate
