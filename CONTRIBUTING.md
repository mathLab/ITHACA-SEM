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
- [Testing and Buildbot](#testing-and-buildbot)
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
6. Respond to any comments in the code review.

## Submission checklist
- Did you add regression tests (for fixes) or unit tests and/or normal tests for
  new features?
- Have you run your branch through buildbot and do all the tests pass?
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

## Testing and Buildbot
Your new features or fixes should include tests that cover the code you've
added. There are numerous examples within the various `Tests` directory lying
within the source trees, and there is an example of writing `.tst` files for our
`Tester` executable in the `tests/Examples` directory. Once you've written your
tests, add them to the `CMakeLists.txt` file for the relevant solver, or to the
appropriate demos directory for library features in whatever directory you are
working in.

You should also test your branch on the
[Nektar++ buildbot](http://buildbot.nektar.info/), which will compile and test
the code against a number of Linux, Mac and Windows operating systems, both 32-
and 64-bit. If your tests don't pass, we can't merge the code into master.

Testing is presently manually executed. You should:

1. Go to the buildbot site and navigate to the *Builders* page.
2. Scroll to the bottom of the page in the section *Force all builds*
3. Enter your details. If you're working on a fork, then the *Suffix to repo
   url* box should be changed to `username/nektar`.
4. Hit the *Force build* button.
5. Check the output in the *Grid* page -- hopefully everything should be green!
   Tests can take up to two hours to run.

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
