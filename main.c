#include <stdio.h>
#include <orbit.h>

extern SatElem OE;
extern double rr[3], vv[3];

int main(int argc, char *argv[])
{
	printf("Hello World!\n");

	rr[0] = 6344.59911  /6378.135;
	rr[1] = -2614.65233 /6378.135;
	rr[2] = -0.00821	/6378.135;

	vv[0] = -0.37178 * 60 / 6378.135;
	vv[1] = -0.88729 * 60 / 6378.135;
	vv[2] = 7.55804	 * 60 / 6378.135;

	rv2el(rr, vv);

	double inc		= OE.xincl,
			xmo		= OE.xmo,
			xno		= OE.xno,
			xnodeo	= OE.xnodeo,
			omegao	= OE.omegao,
			eo		= OE.eo;

	printf("inc: %g		\n", inc / de2ra);
	printf("xmo: %g		\n", xmo / de2ra);
	printf("xno: %g		\n", xno / nocon);
	printf("xnodeo: %g	\n", xnodeo / de2ra);
	printf("omegao: %g	\n", omegao / de2ra);
	printf("eo: %g		\n", eo);

	return 0;
}
