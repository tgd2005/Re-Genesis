ROOT="/Users/murphy92/Y3_Group_Project/oxDNA/oxview-assemblies"

find "$ROOT" -type l | while read -r link; do
  target=$(readlink "$link")

  case "$target" in
    */links/*)
      suffix="${target#*/links/}"
      new_abs="$ROOT/links/$suffix"
      new_rel=$(python3 -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.path.dirname(sys.argv[2])))' "$new_abs" "$link")

      ln -snf "$new_rel" "$link"
      printf 'fixed: %s -> %s\n' "$link" "$new_rel"
      ;;
  esac
done
