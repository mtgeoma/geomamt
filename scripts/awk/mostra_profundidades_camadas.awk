# uso: awk -f mostra_profundidades_camadas.awk z.dat
# z.dat: 1o. coluna=espessura da grade
# z.dat: 2o. coluna=numero da camada a que a grade pertence
{
    h2+=$1;
    printf "%2d [%2d]: %6d %6d\n",NR,$2,h1,h2;h1+=$1;
}
