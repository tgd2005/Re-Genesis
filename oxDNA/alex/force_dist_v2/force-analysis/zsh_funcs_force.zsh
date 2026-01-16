fvec_mag () {
  awk '
    $0 ~ /#/ && $1 ~ /^[0-9]+$/ { t=$1; need_vec=1; next }
    need_vec && NF==3 {
      mag=sqrt($1*$1+$2*$2+$3*$3)
      print t, mag
      need_vec=0
    }
  ' "$@"
}
