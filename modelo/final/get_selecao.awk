# uso: awk -f get_selecao.awk bor004a/xy_bor004a_NM > bor004a/selecao.xy
{
  if(NF==7) {
    printf"%s %s\n",$6,$7
  }
}
