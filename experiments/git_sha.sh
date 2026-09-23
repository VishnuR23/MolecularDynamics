# Shared revision stamp for every artefact this project commits.
#
# The stamp answers one question: which revision of the CODE produced this
# file? So `git describe --always` for the commit, plus a "-dirty" suffix
# only when tracked files outside results/ differ.
#
# Excluding results/ is the whole point. These experiments write into the
# repository, so regenerating one dirties the tree for every experiment that
# runs after it -- and without the exclusion, every run but the first would
# be stamped dirty for a reason that has nothing to do with the code. The
# exclusion keeps "-dirty" meaning what a reader assumes it means: the
# engine or the scripts were modified and this result may not correspond to
# any committed revision.
moldyn_git_sha() {
    local root="$1"
    local base
    base="$(git -C "$root" describe --always 2>/dev/null)" || { printf 'unknown\n'; return; }
    [ -n "$base" ] || { printf 'unknown\n'; return; }
    if [ -n "$(git -C "$root" status --porcelain -- . ':(exclude)results' 2>/dev/null)" ]; then
        printf '%s-dirty\n' "$base"
    else
        printf '%s\n' "$base"
    fi
}
