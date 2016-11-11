/*Copyleft Ashok Palaniappan (ashok@svce.ac.in) and Eric Jakobsson (jake@illinois.edu). GNU GPLv3 Licence.*/
/*Program to calculate moments of a univariate function, here the conservation values of columns of an alignment. */
/*periodicity 1.5 to 6. n = 1 to input value from 5 to 31. output moments at all periodicities, and identify the periodicity with maximum moment.read in values of numerical conservation scores for each position. */

#define MAXLEN_OF_ELEMENT 33
#define MAXLEN_OF_SEQUENCE 1000
#define MAXLEN2_OF_SEQUENCE 100
#define MAXWIDTH_OF_SCORE 16
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/* total #periodicities = (max_period - min_period)*10 +1  = (5.0-2.0)*10+1 =31 */

main (int argc, char *argv[])
{ FILE *infile, *outfile;
int l, m, n, i, temp, eofindicator, seglen, z, cindex, vindex, maxconsmom2pos, maxvarmom2pos, maxconsmom2period, maxvarmom2period, cstartpos[31*MAXLEN2_OF_SEQUENCE], vstartpos[31*MAXLEN2_OF_SEQUENCE];
double period, res_score[MAXLEN_OF_ELEMENT], consmom[31][MAXLEN_OF_SEQUENCE], maxconsmom, maxvarmom, maxconsmom2, maxvarmom2, cperiodicity[31*MAXLEN2_OF_SEQUENCE], vperiodicity[31*MAXLEN2_OF_SEQUENCE], x_comp[31][MAXLEN_OF_SEQUENCE], y_comp[31][MAXLEN_OF_SEQUENCE], maxconsmom2periodicity, maxvarmom2periodicity;
char s[MAXWIDTH_OF_SCORE], c;
const double pi=3.14159;

if (argc!=4) {
printf("Calculation of Conservation Moments\nError: Usage - consmomcalc length_of_segment <input with conservation scores> <outfile>\n");
exit(1);
}

infile=fopen(argv[2],"r");
outfile=fopen(argv[3],"w");

fflush(infile);
rewind(infile);

i=0;
while (1) {
	while (1){
		fscanf(infile,"%s", s);
		if (feof(infile)) {
			eofindicator=1;
			break;
		}
		c=s[0];
		if (!isdigit(c)) continue;
		else break;
	}
if (eofindicator==1) break; 
res_score[i]=atof(s);
i++;
}

seglen=atoi(argv[1]);
/*changing max period; init: (float)(seglen/2); now: 5.0 */

for (period=2.0, l=0; period<=5.0; period+=0.1,l++){
	for (m=0; m<i-seglen+1; m++){
		for (n=0; n<seglen;n++) {

			x_comp[l][m]=x_comp[l][m] + res_score[n+m]*sin(2*pi/period*(n+1));
			y_comp[l][m]=y_comp[l][m] + res_score[n+m]*cos(2*pi/period*(n+1));

		}
	consmom[l][m] = sqrt(pow(x_comp[l][m],2)+pow(y_comp[l][m],2));
	}
}

maxconsmom=0;
maxvarmom=331;
fprintf(outfile, "Variation moments at various periodicities and respective positions for segments of length %d\n\n", seglen);

for (temp=0; temp<l;temp++){
fprintf(outfile, "Periodicity: \t%.1f\n", (float) 2.0 + 0.1*(temp));
	for(m=0;m<i-seglen+1;m++){
		if (consmom[temp][m]>maxconsmom) maxconsmom=consmom[temp][m];
		if (maxvarmom>consmom[temp][m]) maxvarmom=consmom[temp][m];
		fprintf(outfile, "%d. %f\t", m+1, consmom[temp][m]);
	}
	fprintf(outfile, "\n");
}
maxconsmom2=0;
for (temp=0; temp<l;temp++){
	for(m=0;m<i-seglen+1;m++){
		if (consmom[temp][m]==maxconsmom) continue;
		if (consmom[temp][m]>maxconsmom2) {
			maxconsmom2=consmom[temp][m];
			maxconsmom2period=temp;
			maxconsmom2pos=m+1;
		}
	}
}
maxconsmom2periodicity=(double) 2.0+0.1*maxconsmom2period;

maxvarmom2=331;
for (temp=0; temp<l;temp++){
	for(m=0;m<i-seglen+1;m++){
		if (consmom[temp][m]==maxvarmom) continue;
		if (maxvarmom2>consmom[temp][m]) {
			maxvarmom2=consmom[temp][m];
			maxvarmom2period=temp;
			maxvarmom2pos=m+1;
		}
	}
}
maxvarmom2periodicity=(double)2.0+0.1*maxvarmom2period;

cindex=0;
vindex=0;
for (temp=0; temp<l;temp++){
	for(m=0;m<i-seglen+1;m++){
		if (consmom[temp][m]==maxconsmom) {
			cstartpos[cindex]=m+1;
			cperiodicity[cindex]= (double) 2.0 + 0.1*(temp);
			cindex++;
		}
		if (consmom[temp][m]==maxvarmom) {
			vstartpos[vindex]=m+1;
			vperiodicity[vindex] = (double) 2.0+0.1*(temp);
			vindex++;
		}
	}
}

fprintf(outfile, "\n\n");

/* Maximum conservation moment information */
fprintf(outfile, "Maximum Conservation moment value is \t%f.\nIt occurs in %d co-ordinates of (starting position, periodicity).\nCorresponding starting positions, number of residues per turn, and the angular spread per residue in radian and degrees are:\n", maxconsmom, cindex);
for (z=0; z<cindex; z++) {
	fprintf(outfile, "%d\t%f\t%f\t%f\n", cstartpos[z], cperiodicity[z], (float) 2*M_PI/cperiodicity[z], (float) 360/cperiodicity[z]);
}

/*2nd maximum conservation moment*/
fprintf(outfile, "The second maximum conservation moment value is \t%f.\nIt last occurs at position %d along the sequence, with periodicities,  %f residues per turn, %f radians per residue, and %f degrees per residue.\n", maxconsmom2, maxconsmom2pos, maxconsmom2periodicity, (float) 2*M_PI/maxconsmom2periodicity, (float) 360/maxconsmom2periodicity);

/* Maximum Variability moment information */
fprintf(outfile, "Maximum Variability moment value is \t%f.\nIt occurs in %d co-ordinates of (starting position, periodicity).\nCorresponding starting positions, number of residues per turn, and the angular spread per residue in radian and degrees are:\n", maxvarmom, vindex);
for (z=0; z<vindex; z++) {
	fprintf(outfile, "%d\t%f\t%f\t%f\n", vstartpos[z], vperiodicity[z], (float) 2*M_PI/vperiodicity[z], (float) 360/vperiodicity[z]);
}

/*2nd maximum Variability moment*/
fprintf(outfile, "The second maximum variability moment value is \t%f.\nIt last occurs at position %d along the sequence, with periodicities,  %f residues per turn, %f radians per residue, and %f degrees per residue.\n", maxvarmom2, maxvarmom2pos, maxvarmom2periodicity, (float) 2*M_PI/maxvarmom2periodicity, (float) 360/maxvarmom2periodicity);

/* calculation of strength of helical and beta-shhet fourier components*/

fprintf(outfile,"\nConservation at helical periodicities\n");
for(m=0;m<i-seglen+1;m++){
	fprintf(outfile,"%d\t%f\n", m+1, consmom[16][m]);
}

fprintf(outfile, "\nConservation with beta-structures\n");
for(m=0;m<i-seglen+1;m++){
	fprintf(outfile,"%d\t%f\n", m+1, consmom[0][m]);
}


}
