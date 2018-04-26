#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <orbit.h>
#include "csvparser.h"

extern SatElem OE;

///////////// declare global variables ////////////////////////////////////////

double tle, inc, xnodeo, eo, omegao, xmo, xno;
char buf[21], name[20],
	 sn[] = "12345", dsg[] = "12345A", bstar[] = "10000-4";
extern double rr[3], vv[3];

// write TLE to output file and to screen
void write_tle(char *file_out)
{
  char eo_string[7];
  char line1[70];
  char line2[70];

  sprintf(eo_string, "%.7f", eo);
  eo_string[0] = eo_string[2];
  eo_string[1] = eo_string[3];
  eo_string[2] = eo_string[4];
  eo_string[3] = eo_string[5];
  eo_string[4] = eo_string[6];
  eo_string[5] = eo_string[7];
  eo_string[6] = eo_string[8];

  double xndt2o = 0;

  sprintf(line1, "1 %-.5sU %-.6s   0%14.8f %.8f  00000-0   %7s 0    00"
		   ,sn, dsg, tle, xndt2o, bstar);
   // change tle format character from 0%13.8f to %14.8f on 1/1/2010
//  ccksum(line1);
  sprintf(line2, "2 %-.5s %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
	sn, inc, xnodeo, eo_string, omegao, xmo, xno);
//  ccksum(line2);
  FILE *fp;
  fp = fopen(file_out, "VEC_KM.txt");
  fprintf(fp, "\n%s\n", name);
  fprintf(fp, "%s\n", line1);
  fprintf(fp, "%s\n", line2);
  fclose(fp);
  printf("\n%s\n", name);
  printf("%s\n", line1);
  printf("%s\n", line2);
}

int main(int argc, char *argv[])
{
    int i =  0;
    //                                   file, delimiter, first_line_is_header?
	CsvParser *csvparser = CsvParser_new("table.csv", ";", 1);
    CsvRow *header;
    CsvRow *row;

    header = CsvParser_getHeader(csvparser);
    if (header == NULL) {
        printf("%s\n", CsvParser_getErrorMessage(csvparser));
        return 1;
    }
    char **headerFields = CsvParser_getFields(header);
    for (i = 0 ; i < CsvParser_getNumFields(header) ; i++) {
        printf("TITLE: %s\n", headerFields[i]);
    }
    // CsvParser_destroy_row(header); -> causes error in current version
    while ((row = CsvParser_getRow(csvparser)) ) {
        printf("NEW LINE:\n");
        char **rowFields = CsvParser_getFields(row);
        for (i = 0 ; i < CsvParser_getNumFields(row) ; i++) {
            printf("FIELD: %s\n", rowFields[i]);
			if (i < 3)
				rr[i] = atof(rowFields[i]) / 6378.135;
			if (i >= 3)
				vv[i - 3] = atof(rowFields[i]) * 60 / 6378.135;
        }
		for (i = 0 ; i < 3 ; i++){
			printf("%8.4lf\n", rr[i]);
		}
		for (i = 0 ; i < 3 ; i++){
			printf("%8.4lf\n", vv[i]);
		}
		rv2el(rr, vv);

		inc		= OE.xincl;
		xmo		= OE.xmo;
		xno		= OE.xno;
		xnodeo	= OE.xnodeo;
		omegao	= OE.omegao;
		eo		= OE.eo;

		printf("\ninc: %g		\n", inc / de2ra);
		printf("xmo: %g		\n", xmo / de2ra);
		printf("xno: %g		\n", xno / nocon);
		printf("xnodeo: %g	\n", xnodeo / de2ra);
		printf("omegao: %g	\n", omegao / de2ra);
		printf("eo: %g		\n", eo);

        CsvParser_destroy_row(row);
    }
    CsvParser_destroy(csvparser);



//	char file_in[20], file_out[20];

//	sprintf(file_in, "VEC_KM.txt");
//	sprintf(file_out, "VEC_KM.txt");

//	for (int i = 0; i < 3; i++) {
//		printf("\nEnter rr[%d]: ", i);
//		scanf("%lf", &rr[i]);
//	}

//	for (int i = 0; i < 3; i++) {
//		printf("\nEnter vv[%d]: ", i);
//		scanf("%lf", &vv[i]);
//	}

//	for (int i = 0; i < 3; i++) {
//		printf("\nrr[%d]= %lf\nvv[%d]= %lf", i, rr[i], i, vv[i]);
//	}

//	rr[0] /= 6378.135;		// 4625.95913884436
//	rr[1] /= 6378.135;		// -280.245682707435
//	rr[2] /= 6378.135;		// 5043.35624112831

//	vv[0] = vv[0] * 60 / 6378.135;	// -5.59865049836795
//	vv[1] = vv[1] * 60 / 6378.135;	// -1.08162703917961
//	vv[2] = vv[2] * 60 / 6378.135;	// 5.06186754556158

//	for (int i = 0; i < 3; i++) {
//		printf("\nrr[%d]= %lf\nvv[%d]= %lf", i, rr[i], i, vv[i]);
//	}

//	rv2el(rr, vv);

//	inc		= OE.xincl;
//	xmo		= OE.xmo;
//	xno		= OE.xno;
//	xnodeo	= OE.xnodeo;
//	omegao	= OE.omegao;
//	eo		= OE.eo;

////	printf("\ninc: %g		\n", inc / de2ra);
////	printf("xmo: %g		\n", xmo / de2ra);
////	printf("xno: %g		\n", xno / nocon);
////	printf("xnodeo: %g	\n", xnodeo / de2ra);
////	printf("omegao: %g	\n", omegao / de2ra);
////	printf("eo: %g		\n", eo);

//	write_tle(file_out);

	return 0;
}
