# Releasing ymf

A release is a git tag. `.github/workflows/release.yml` builds the sdist and
wheel from the tagged commit, checks them (`tools/check_dist.sh`: size,
contents, `twine check`, a fresh-venv install that runs), and publishes with
PyPI Trusted Publishing (OIDC: no API token, no secret). The version is
whatever `setuptools_scm` derives from the tag; the workflow refuses a tag that
does not equal `v<built version>`.

| tag | goes to |
|---|---|
| `v0.3.0rc1`, `v0.3.0a1`, `v0.3.0b1` | TestPyPI (environment `testpypi`) |
| `v0.3.0` | PyPI (environment `pypi`) |

The tag `v0.2.0` already exists from before this workflow; it is not re-pushed
and nothing publishes it. Pick the first published version knowing that.

## What the maintainer does by hand (once)

Name check, 2026-11-16: `ymf` returns 404 on pypi.org and test.pypi.org (JSON
API and simple index), so it is free today. PyPI can still refuse a name that
is confusingly close to an existing project, and a project registered
meanwhile wins: do step 1 soon.

1. **TestPyPI**: log in at test.pypi.org (a separate account from PyPI), then
   *Your account → Publishing → Add a new pending publisher*: project name
   `ymf`, owner `proteustoolkit` (or `cekees` to rehearse from the fork),
   repository `ymf`, workflow `release.yml`, environment `testpypi`.
2. **PyPI**: the same at pypi.org, environment `pypi`. A *pending* publisher
   creates the project on the first upload, so nothing is uploaded to claim
   the name.
3. **GitHub environments** on the repository that will run the workflow:
   Settings → Environments → create `testpypi` and `pypi`. For `pypi`, add
   yourself as a required reviewer so a final tag waits for a click; leave
   `testpypi` unprotected. No secrets are added to either.
4. **Push a pre-release tag** from a commit on `main` whose CI is green:

   ```bash
   git tag -a v0.3.0rc1 -m "ymf 0.3.0rc1" && git push origin v0.3.0rc1
   ```

   Then try it the way a user would:

   ```bash
   python -m venv /tmp/t && /tmp/t/bin/pip install --index-url https://test.pypi.org/simple/ \
       --extra-index-url https://pypi.org/simple/ "ymf[frontend]==0.3.0rc1"
   ```

   (the extra index supplies pyyaml and the extras' dependencies, which are
   not on TestPyPI). A TestPyPI file name can never be reused: a mistake is
   fixed by `rc2`, not by re-uploading.
5. When rc installs and runs, tag `v0.3.0` the same way for PyPI.

To rehearse the build without publishing, run `tools/check_dist.sh` (set
`SETUPTOOLS_SCM_PRETEND_VERSION=0.3.0rc1` off a tag).

## Pinning from proteus CI (not done here)

proteus's workflows install ymf from `git+https://github.com/proteustoolkit/ymf.git@main`,
a moving branch, in the environment files and in the comment block of
`.github/workflows/sdist.yml`. After a tagged release exists, the proteus-side
change is: replace `@main` with the tag (`@v0.3.0`, or `ymf==0.3.0` from PyPI
once it is there; before that, TestPyPI needs the extra-index form above and
is not worth using in CI), in every environment file and workflow that names
it, and update the comment. The pin moves by a one-line PR per ymf release; a
ymf change that proteus needs merges in ymf, is tagged, and then the proteus PR
bumps the pin. That is a proteus change: it needs its own PR and is not part of
this one.
