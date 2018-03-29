#include <stdio.h>
#include <stdlib.h>                 /*  Canny */
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

typedef enum {false, true} bool;
int    pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double Histogram[PICSIZE] = { 0 };
double mag[PICSIZE][PICSIZE];
double peak[PICSIZE][PICSIZE];
double fpeaks[PICSIZE][PICSIZE];
bool peaks[PICSIZE][PICSIZE];
bool flags[PICSIZE][PICSIZE];

int main(argc,argv)
int argc;
char **argv;
{
    int     i,j,p,q,mr,centx,centy, sum1, sum2, sig;
    double  xmaskval, ymaskval, maxival, slope, percent, CutOff, AreaOfTops=0, HI, LO;
    FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
    char    *foobar;

    // input image
    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");
    // magnitude output image
    argc--; argv++;
    foobar = *argv;
    fo1=fopen(foobar,"wb");
    // peaks output image
    argc--; argv++;
    foobar = *argv;
    fo2=fopen(foobar,"wb");
    // final algorithm output image
    argc--; argv++;
    foobar = *argv;
    fo3=fopen(foobar,"wb");
    // sigma value
    argc--; argv++;
    foobar = *argv;
    sig = atoi(foobar);
    // percent value
    argc--; argv++;
    foobar = *argv;
    percent = atof(foobar);

    // formatting for .pgm
    fprintf(fo1,"P5\n");
    fprintf(fo1,"%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo1,"255\n");
    fprintf(fo2,"P5\n");
    fprintf(fo2,"%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo2,"255\n");
    fprintf(fo3,"P5\n");
    fprintf(fo3,"%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo3,"255\n");

    // setting the mask radius and centers of x and y
    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    // translating input image to array
    for (i=0;i<PICSIZE;i++)
    { for (j=0;j<PICSIZE;j++)
            {
              pic[i][j]  =  getc (fp1);
            }
    }

    // getting the mask for x and y through partial derivatives
    for (p=-mr;p<=mr;p++)
    {  for (q=-mr;q<=mr;q++)
       {
          xmaskval = ((-q)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
          ymaskval = ((-p)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
          maskx[p+centy][q+centx] = xmaskval;
          masky[p+centy][q+centx] = ymaskval;
       }
    }

    // convolution
    for (i=mr;i<PICSIZE-mr;i++)
    { for (j=mr;j<PICSIZE-mr;j++)
      {
         sum1 = 0;
         sum2 = 0;
         for (p=-mr;p<=mr;p++)
         {
            for (q=-mr;q<=mr;q++)
            {
               sum1 += pic[i+p][j+q] * maskx[p+centy][q+centx];
               sum2 += pic[i+p][j+q] * masky[p+centy][q+centx];
            }
         }
         outpicx[i][j] = sum1;
         outpicy[i][j] = sum2;
      }
    }

    // getting the magnitude array through root of squares
    maxival = 0;
    for (i=mr;i<PICSIZE-mr;i++)
    { for (j=mr;j<PICSIZE-mr;j++)
      {
         mag[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                  (outpicy[i][j]*outpicy[i][j])));
         if (mag[i][j] > maxival)
            maxival = mag[i][j];

       }
    }

    // print magnitude image
    for (i=0;i<PICSIZE;i++)
    {
        for (j=0;j<PICSIZE;j++)
        {
         mag[i][j] = (mag[i][j] / maxival) * 255;
         fprintf(fo1,"%c",(char)((int)(mag[i][j])));

        }
    }

    // calculating which pixels are peaks with the canny algorithm
    for (i=mr;i<PICSIZE-mr;i++)
    {
        for (j=mr;j<PICSIZE-mr;j++)
        {
            if(outpicx[i][j] == 0.0)
            {
                outpicx[i][j] = .00001;
            }
            slope = outpicy[i][j]/outpicx[i][j];
            if((slope <= .4142) && (slope > -.4142))
            {
                if((mag[i][j] > mag[i][j-1]) && (mag[i][j] > mag[i][j+1])){
                    peak[i][j] = 255.0;
                }
            }
            else if((slope <= 2.4142) && (slope > .4142))
            {
                if((mag[i][j] > mag[i-1][j-1]) && (mag[i][j] > mag[i+1][j+1])){
                    peak[i][j] = 255.0;
                }
            }
            else if((slope <= -.4142) && (slope > -2.4142))
            {
                if((mag[i][j] > mag[i+1][j-1]) && (mag[i][j] > mag[i-1][j+1])){
                    peak[i][j] = 255.0;
                }
            }
            else
            {
                if((mag[i][j] > mag[i-1][j]) && (mag[i][j] > mag[i+1][j])){
                    peak[i][j] = 255.0;
                }
            }
        }
    }

    // printing new peaks array
    for (i=0;i<PICSIZE;i++)
    {
        for (j=0;j<PICSIZE;j++)
        {
         fprintf(fo2,"%c",(char)((int)(peak[i][j])));

        }
    }

    // making a histogram of the magnitude array
    for (i=0;i<PICSIZE;i++)
    {
        for (j=0;j<PICSIZE;j++)
        {
            Histogram[(int)mag[i][j]]++;
        }
    }

    CutOff = percent*PICSIZE*PICSIZE;

    for(i=255; i>=0;i--)
    {
        AreaOfTops += Histogram[i];
        if(AreaOfTops > CutOff)
            break;

    }

    // thresholds
    HI=i;
    LO=(.35*HI);

    // initializing arrays for iterative sweep of the peaks array
    for (i=0;i<PICSIZE;i++)
    {
        for (j=0;j<PICSIZE;j++)
        {
            if(peak[i][j] == 255.0)
            {
                // if the peak is above the high threshold, automatically goes
                // to the final output image
                peaks[i][j] = true;
                if(mag[i][j] > HI)
                {
                    peaks[i][j] = false;
                    flags[i][j] = true;
                }
                else if(mag[i][j] < LO)
                {
                    peaks[i][j] = flags[i][j] = false;
                }
            }
        }
    }

    // go through peak array to find which pixels are in between the high and
    // low thresholds to determine if they will go to the final output image
    // by relating them to their neighbors
    bool moretodo = true;
    while(moretodo == true)
    {
        moretodo =false;
        for (i=0;i<PICSIZE;i++)
        {
            for (j=0;j<PICSIZE;j++)
            {
                if(peaks[i][j] == true)
                {
                    for(p=-1; p<=1; p++)
                    {
                        for(q=-1; q<=1; q++)
                        {
                            if(flags[i+p][j+q] == true)
                            {
                                moretodo = true;
                                flags[i][j] = true;
                                peaks[i][j] = false;
                            }
                        }
                    }
                }
            }
        }
    }

    // print new final canny peak array
    for (i=0;i<PICSIZE;i++)
    {
        for (j=0;j<PICSIZE;j++)
        {
         if(flags[i][j] == true) fpeaks[i][j] = 255;
         fprintf(fo3,"%c",(char)((int)(fpeaks[i][j])));

        }
    }

    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    fclose(fp1);
}


