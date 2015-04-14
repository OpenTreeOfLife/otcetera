Subproblems are written in the export-sub-temp directory here, and 
then moved to their final location only if they differ from content there
(based on diffing their checksums).

This prevents touching subproblems that have not changed, so that 
solutions for unchanged subproblems are not needlessly recalculated.