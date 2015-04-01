# usos:
# echo "<rho> <T>" | awk -f skin-depth.awk
# espera que a 1a. coluna seja a resistividade aparente [ohm.m] e a 2a. o período [s]
# saída= skin-depth [m]
#
# echo "<sd> <rho>" | awk -f skin-depth.awk -v d=T
# espera que a 1a. coluna seja o skin-depth [m] e a 2a. a resistividade aparente [ohm.m]
# saída= T [s]
#
# echo "<sd> <T>" | awk -f skin-depth.awk -v d=rho
# espera que a 1a. coluna seja o skin-depth [m] e a 2a. o período [s]
# saída= resistividade aparente [ohm.m]
#
# sd=sqrt((rho*T)/(pi*u));
BEGIN{
    pi=atan2(1,1)*4;
    u=4e-7*pi;
}
{
    if(d=="T"){
	sd=$1;
	rho=$2;
	printf "%f\n",(sd*sd)*(pi*u)/rho;
    }
    else if(d=="rho") {
	sd=$1;
	T=$2;
	printf "%f\n",(sd*sd)*(pi*u)/T;
    }
    else {
	rho=$1;
	T=$2;
	printf "%f\n",sqrt((rho*T)/(pi*u));
    }
}
