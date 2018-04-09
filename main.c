#include <stdio.h>
#include <orbit.h>

extern OE_t OE;
double rr[3], vv[3];

int main(int argc, char *argv[])
{
	printf("Hello World!\n");

	rr[0] = 6344.59911;
	rr[1] = -2614.65233;
	rr[2] = -0.00821;

	vv[0] = -0.37178;
	vv[1] = -0.88729;
	vv[2] = 7.55804;

	rv2el(rr, vv);

	printf("struct: ", OE.xincl);

	return 0;
}
