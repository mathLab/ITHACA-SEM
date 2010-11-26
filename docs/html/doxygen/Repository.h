/**
 * \page pageRepositoryUsage SVN Repository Usage
 * Here we detail guidance on using the subversion repository. In particular,
 * we describe how to use branches for feature development and how to commit
 * the branches back into the trunk when completed.
 *
 * \section sectionRepositoryUsageStructure Repository structure
 * The root of the Nektar++ repository is
 * https://gforge.sci.utah.edu/svn/nektar. The standard convention for code
 * repositories stipulates three sub-directories: trunk, tags, branches. The
 * trunk directory is the primary version of the code which should contain the
 * latest completed features. The branches directory contains (temporary)
 * copies of trunk in which new features can be developed without introducing
 * instability into the trunk. Once the feature is complete and merged back into
 * the trunk, the branch is usually deleted. The tags directory contains
 * snapshots of the trunk at notable points in time (e.g. release candidates).
 *
 * Revision numbers are repository-global and incremented whenever a commit is
 * made in a trunk, a branch or in creating a tagged version. Therefore,
 * the sequence of revision numbers for a given branch or the trunk will not
 * necessarily be consecutive.
 *
 * \section sectionRepositoryUsageBranches Branches
 * Branches are complete copies of the main (trunk) code tree and allow
 * parallel development of partially complete features without them interfering
 * with each other. When the feature is complete, the changes are merged back
 * into the trunk. In this way, the trunk remains stable.
 *
 * \subsection sectionRepositoryUsageBranchesCreate Creating a branch
 * To create a branch:
 * @code
 * svn copy https://gforge.sci.utah.edu/svn/nektar/trunk \
 *   https://gforge.sci.utah.edu/svn/nektar/branches/my-feature
 * @endcode
 * where my-feature is a name for the branch.
 *
 * \subsection sectionRepositoryUsageBranchesMerge Merge a branch back into trunk
 * We will assume we have a branch 'my-feature' which we wish to merge into the
 * trunk. Before starting, determine the following items of information:
 * - Revision at which the branch was created (denote as 'M'). From a working
 *   copy of your branch run:
 *   @code
 *   svn log --stop-on-copy
 *   @endcode
 *   and look for the most recent entry in which code from the trunk is copied
 *   into the branch (e.g. branch creation or the latest merge).
 * - The latest revision number of the repository. From either trunk or branch
 *   @code
 *   svn update
 *   svn info
 *   @endcode
 *   and look at the 'Revision' field. Denote this as 'N'.
 *
 * We merge the branch into the trunk and (optionally) pull any changes to
 * the trunk since the branch was created into the branch. The
 * latter is only necessary if we plan to continue working in the branch.
 * @note In this case, for subsequent merges we will pull in changes since our
 * last merge, rather than the branch creation.
 * @endnote
 *
 * We begin by merging the changes to our branch into the trunk.
 * -# Ensure everything in the branch is committed. From your working copy of
 *    the branch run
 *    @code
 *    svn status
 *    @endcode
 *    If there are any files listed with anything other than a ? next to them,
 *    you should commit these changes if necessary.
 * -# Check out a fresh working copy of the latest trunk
 *    @code
 *    svn checkout https://gforge.sci.utah.edu/svn/nektar/trunk nek-trunk
 *    @endcode
 * -# From the working copy of the trunk, perform a dry-run of the merge:
 *    @code
 *    svn merge --dry-run -r M:HEAD https://gforge.sci.utah.edu/svn/nektar/branches/my-feature
 *    @endcode
 *    replacing M with the start revision you noted down earlier. All the
 *    changed files should be listed. If any have a C next to them there is a
 *    conflict between the changes you made in your branch and other changes
 *    which have been made in the trunk. If there are a lot of C's, it is
 *    likely you have done something wrong with the instructions above.
 * -# If everything looks good, perform the merge:
 *    @code
 *    svn merge -r M:HEAD https://gforge.sci.utah.edu/svn/nektar/branches/my-feature
 *    @endcode
 * -# Rebuild everything in the trunk, run regression tests, check everything
 *    works as expected. Test your new feature works, and that existing features
 *    in the trunk still work.
 * -# Commit the trunk to apply the merged changes to the repository.
 *    @code
 *    svn update
 *    svn commit
 *    @endcode
 *    It is important to give details of exactly what was merged in this commit.
 *    Always give a commit message such as
 *    @code
 *    "Merged revisions M:N of branch my-feature into trunk"
 *    @endcode
 *    replacing 'M', 'N' and 'my-feature' as appropriate.
 *
 * Optionally, we may now incorporate changes to the trunk which have happened
 * since the branch was created into our branch to continue using it.
 * -# Change to your working copy of your branch.
 * -# Ensure it is up-to-date
 *    @code
 *    svn update
 *    @endcode
 * -# Merge changes between branch creation and revision N (i.e. the revision
 *    number of the merge commit above, minus one). First do a dry-run.
 *    @code
 *    svn merge --dry-run -r M:N https://gforge.sci.utah.edu/svn/nektar/trunk
 *    @endcode
 *    Check the proposed changes look correct. Many files with C might indicate
 *    an incorrect choice of M and N.
 * -# Perform the actual merge:
 *    @code
 *    svn merge -r M:N https://gforge.sci.utah.edu/svn/nektar/trunk
 *    @endcode
 * -# Rebuild your branch, run regression tests, check everything works.
 * -# Commit the modifications to your branch.
 *    @code
 *    svn commit
 *    @endcode
 *    Again, always give a commit message, explaining exactly what was merged.
 *
 * As a final check, the trunk and your branch should now be identical. This
 * can be checked using:
 * @code
 * svn diff https://gforge.sci.utah.edu/svn/nektar/trunk \
 *   https://gforge.sci.utah.edu/svn/nektar/branches/my-feature
 * @endcode
 */
