Contributing to Nektar++
========================

## Contents
This is a reasonably complete guide to help if you're interested in contributing
to Nektar++, either in reporting bugs or, hopefully, trying to fix them! It's
split up into a number of sections:

- [Issues and bug reports](#issues-and-bug-reports)
- [How to contribute](#how-to-contribute)
- [Submission checklist](#submission-checklist)
- [Git cheatsheet](#git-cheatsheet)
- [Testing and GitLab CI](#testing-and-gitlab)
- [Documentation](#documentation)
- [Formatting guidelines](#formatting-guidelines)

## Issues and bug reports
Think you've found a bug or issue with Nektar++? We're very keen to hear about
it!
- In the first instance, you should raise an issue on the
  **[issue tracker](https://gitlab.nektar.info/nektar/nektar/issues)** -- be
  sure to do a quick search and see if anyone has reported the same thing first.
- Alternatively you can
  **[join the mailing list](https://mailman.ic.ac.uk/mailman/listinfo/nektar-users)**
  for more advice.

It's *really helpful* if you can include a small session file that reproduces
the error, and can give a good description of the problem you're having.

## How to contribute
If you've got a patch or feature, please consider contributing it back to the
project. It's a pretty simple process:

1. Fork the Nektar++ repository in `nektar/nektar` into your username's space.
2. Create a branch with the naming convention:
   - `feature/myawesomebranch`: a new feature that wasn't in Nektar++ already.
   - `fix/mygreatfix`: fixes an issue that isn't tracked in the issue tracker.
   - `ticket/123-myfantasticpatch`: fixes an issue that is tracked in the issue
     tracker (please include the issue number somewhere!)
   - `tidy/mybrillianttidying`: cosmetic fixes to bring existing files up to the
     Nektar++ code guidelines.
3. Make sure you've gone through the checklist below.
4. Submit a merge request to merge into `master`. If you just want to see the
   diff and are not quite ready to merge, use the `[WIP]` tag in the title to
   prevent your code from being accidentally merged.
5. Put a comment in the MR saying that it's ready to be merged.
6. If your branch is a minor fix that could appear in the next patch release,
   then add the `Proposed patch` label to the merge request.
7. Respond to any comments in the code review.

## Submission checklist
- Did you add regression tests (for fixes) or unit tests and/or normal tests for
  new features?
- Have you run your branch through GitLab CI and do all the tests pass?
- Have you fixed any new compiler warnings your code has introduced into the
  compilation step for all of the Linux CI environments?
  - **unused parameters**: if these are genuinely needed (e.g. virtual functions
    in a derived class, please use `boost::ignore_unused()` to mark as such.
  - **switch case may fall-through**: for switch statements which
    *intentionally* exploit fall-through between cases, mark the end of such
    cases with the comment `/* Falls through. */` to suppress the warning.
  - Avoid `ASSERTL0(false, msg)`; instead use `NEKERROR(ErrorUtil:efatal, msg)`.
  - Ensure variables are initialised with sensible default values.
- Is there documentation in the user guide and/or developer guide?
- Have you added a CHANGELOG entry, including the MR number?
- Are there any massive files you might have added in the commit history? We try
  to keep test files as small as possible. If so you'll need to rebase or
  filter-branch to remove those from the commit history.
- Is the code formatted correctly?
  - **Note:** unfortunately, Nektar++ has pretty inconsistent code formatting at
    the moment. To help in reviewing your submission, new files should be
    formatted according to the guidelines (or use `clang-format` as described
    below) -- otherwise, try to keep formatting consistent with the file you're
    working on.

## Git cheatsheet
Although Gitlab gives a nice interface to view the diff between a branch and
master, for large merges, it can be slow. The following `git` aliases can
provide a quicker alternative. You can use these by inserting them into the
`.gitconfig` file in your home directory, or inside the `nektar++/.git/config`
file.

```
[alias]
branch-name = "!git rev-parse --abbrev-ref HEAD"
diff-nows = diff --color -w
log-branch = log --pretty='%C(green)%h %C(red)%an %C(reset)(%C(blue)%ad%C(reset))%n%s' master..
diff-branch = diff -U5 --minimal --color -w master...
```

This gives you four commands:

- `git branch-name` displays the current branch name
- `git diff-nows` shows a diff of your current commit in colour, without
  whitespace changes.
- `git log-branch` shows a minimised log of all the commits on the current
  branch that are not in `master`.
- `git diff-branch` shows a diff of the current branch against `master`, without
  showing changes from `master` that aren't present in the branch (i.e. `git
  diff master...branch`), without whitespace changes. (This should be roughly
  equivalent to Gitlab's diff).

If you prefer a graphical interface to see the files that have changed in your
commit, you can additionally use the `git gui` command to bring up a simple
interface. `git difftool` can also be used in combination with a GUI diff
viewer, to graphically view the output of `git diff`.

## Testing and GitLab CI
Your new features or fixes should include tests that cover the code you've
added. There are numerous examples within the various `Tests` directory lying
within the source trees, and there is an example of writing `.tst` files for our
`Tester` executable in the `tests/Examples` directory. Once you've written your
tests, add them to the `CMakeLists.txt` file for the relevant solver, or to the
appropriate demos directory for library features in whatever directory you are
working in.

You should also test your branch on the Nektar++ GitLab CI, which will compile
and test the code against a number of Linux, Mac and Windows operating
systems. If your tests don't pass, we can't merge the code into master.

When you submit a merge request testing on GitLab CI will happen automatically,
unless you have marked the merge request as a work-in-progress (WIP: prefix).
Each time you push commits to a non-WIP merge request branch, it will also
trigger a build.

## Documentation
Nektar++ has a fairly comprehensive user guide and a developer guide that is
presently very incomplete. The following are rough guidelines for what you
should provide:

- If you are writing user-exposed features, you should add some documentation to
  the user guide on how to use them.
- Any functions/classes should include Doxygen documentation.
- Generally, code should be well-commented using regular C++ comments to explain
  its function to help in reviewing it.

Nektar++ also has a growing number of tutorials to help introduce users and
developers to the use of the library and the range of application solvers. These
are stored in a separate repository, but are available from the main repository
through a git submodule. To populate the docs/tutorial directory run `git
submodule init` followed by `git submodule update --remote`. The latter command
will ensure you have the latest master branch of the tutorials within your
source tree.

## Code review and merging
All merge requests will be reviewed by one of the senior developers. We try to
stick to the following process:
- Senior developer will be assigned, MR will be assigned a milestone to target a
  release.
  - If the branch is deemed to be minor and passes the checklist above, senior
    developer will handle the request by themselves.
  - Otherwise, senior developer will ask one or more other developers to review
    the code.
- Submission checklist will be checked by the reviewers.
- Where appropriate, reviewers will comment on regions of code that need further
  development and/or improvement.
- In addition to any coding comments/suggestions, reviewers are asked to check
  the branch passes the regression tests and appropriate documentation has been
  added.
- Once feedback received from the branch author (if necessary) and reviewers are
  happy, the branch will be merged.

## Release branches
Nektar++ releases are versioned in the standard form `x.y.z` where `x` is a
major release, `y` a minor release and `z` a patch release:

- major releases are extremely infrequent (on the order of every 2-3 years) and
  denote major changes in functionality and the API;
- minor releases occur around twice per year and contain new features with minor
  API changes;
- patch releases are targeted on roughly a monthly basis and are intended to
  fix minor issues in the code.

The repository contains a number of _release branches_ named `release/x.y` for
each minor release, which are intended to contain **fixes and very minor
changes** from `master` and which form the next patch release. This allows us to
use `master` for the next minor release, whilst still having key fixes in patch
releases.

### Cherry-picking process

Any branches that are marked with the `Proposed patch` label should follow the
following additional steps to cherry pick commits into the `release/x.y` branch.

1. If the branch is on a remote other than `nektar/nektar`, make sure that's
   added to your local repository.
2. On a local terminal, run `git fetch --all` to pull the latest changes. It's
   important for the commands below that you do this _before_ you merge the
   branch into `master`.
3. Merge the branch into master as usual using GitLab.
4. Switch to the appropriate branch with `git checkout release/x.y` and update
   with `git pull`.
5. Now check the list of commits to cherry-pick.

   ```bash
   git log --oneline --no-merges --reverse origin/master..REMOTE/fix/BRANCHNAME
   ```

   where `REMOTE` is the remote on which the branch lives and `BRANCHNAME` is
   the fix branch. If the list is empty, you probably did a `git fetch` after
   you merged the branch into `master`; in this case use `origin/master^`.
6. If you're happy with the list (compare to the MR list on the GitLab MR if
   necessary), cherry-pick the commits with the command:

   ```bash
   git cherry-pick -x $(git rev-list --no-merges --reverse origin/master..REMOTE/fix/BRANCHNAME)
   ```

7. It's likely you'll encounter some conflicts, particularly with the
   `CHANGELOG`. To fix these:
   - `git status` to see what's broken
   - Fix appropriately
   - `git commit -a` to commit your fix
   - `git cherry-pick --continue`
8. If everything becomes horribly broken, `git cherry-pick --abort`.
9. Once you're happy, `git push` to send your changes back to GitLab.

Steps 5 and 6 can be simplified by creating a script
```bash
#!/bin/bash
src=$1

logopts="--oneline --no-merges --reverse"
commits=`git log $logopts master..$1 | cut -f 1 -d " " | xargs`

echo "Will cherry-pick the following commits: $commits"
echo "Press ENTER to continue..."
read

cherryopts="-x --allow-empty --allow-empty-message"
git cherry-pick $cherryopts $commits
```
which accepts the name of the source branch as the sole argument.

## Formatting guidelines
Nektar++ uses C++, a language notorious for being easy to make obtuse and
difficult to follow code. To hopefully alleviate this problem, there are a
number of fairly simple formatting guidelines you should follow. We are
reasonably relaxed about code formatting, but if you can follow the guidelines
below this would be fantastic.

### Basic rules
- All code should be wrapped to 80 characters.
- Indentation should be 4 spaces with **no tabs**. Namespaces should not be
  indented to give more room in the 80 character width.
- Please comment your code with Doxygen and inline comments wherever possible --
  but don't use trailing inline comments to save the 80 character limit!
- All code blocks (even one-line blocks) should use braces, and braces should be
  on new lines; for instance
  
  ```c++
  if (someCondition)
  {
      myAwesomeFunction();
  }
  ```

- **Don't use preprocessor directives and macros unless there is no viable
  alternative.**
- However, please make sure you do have a header guard inside your `.h` files,
  which you should be sure to include in any headers you contribute.
- Use one `.cpp` and `.h` file per C++ class, and try to keep `inline` header
  code to a minimum (unless performance is a major factor).
- Put spaces around binary operators and constants.
- Put spaces after `if`, `while`, etc., but not after function names (see the
  example above).

### Variables and naming
- Please use sensible names and use camelCase as a broad naming convention.
  - Variables should start with a lowercase letter, e.g. `myAwesomeVariable`.
  - Function, `class`, `struct` and `typedef` names should begin with capital
    letters, e.g. `MyAwesomeFunction`.
- Inside classes, member variables should be prefixed with `m_`,
  e.g. `m_myAwesomeVariable`.
  - Global constants used throughout the library should be prefixed with `k`
    (e.g. `kGeometricTolerance`), and enumerations should be prefixed with `e`
    (e.g. `eGeometry`).
- Use all uppercase letters with underscores between words for pre-processor
  definitions and macros.

### Using `clang-format`
Code formatting is reasonably boring, so Nektar++ comes with a `.clang-format`
file to allow for automatic code formatting. As noted above, you can use this
for new files, or cosmetic `tidy/*` branches, but try to stick to existing
formatting elsewhere.

Installing it is straightforward on most package managers. Nektar++ relies on
options that are used in version 3.7 or later.

There are a number of instructions on how to use `clang-format` inside a number
of text editors on the
[CLang website](http://clang.llvm.org/docs/ClangFormat.html). However at a
minimum, you should consider downloading the
[``git-clang-format``](https://llvm.org/svn/llvm-project/cfe/trunk/tools/clang-format/git-clang-format)
script into one of your `$PATH` locations. You can then run the command

    git clang-format

before you do a `git commit`, and `clang-format` will automatically format your
diff according to the guidelines.
