#include <stdio.h>
#include <stdlib.h>                 /*  Marr-Hildreth.c  (or marrh.c) */
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
    double Histogram[256] = { 0 };
    double ival[256][256];
    double mag[256][256];
    double peak[256][256];
    double fpeaks[256][256];
    bool peaks[256][256];
    bool flags[256][256];

int main(argc,argv)
int argc;
char **argv;
{
    int     i,j,p,q,mr,centx,centy, sum1, sum2, sig;
    double  xmaskval, ymaskval, maxival, slope, percent, CutOff, AreaOfTops=0, HI, LO;
    FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
    char    *foobar;
    int cols = 256, rows = 256;

    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    argc--; argv++;
    foobar = *argv;
    fo1=fopen(foobar,"wb");

    argc--; argv++;
    foobar = *argv;
    fo2=fopen(foobar,"wb");

    argc--; argv++;
    foobar = *argv;
    fo3=fopen(foobar,"wb");

    argc--; argv++;
    foobar = *argv;
    sig = atoi(foobar);

    argc--; argv++;
    foobar = *argv;
    percent = atof(foobar);

    fprintf(fo1,"P5\n");
    fprintf(fo1,"%d %d\n", rows, cols);
    fprintf(fo1,"255\n");

    fprintf(fo2,"P5\n");
    fprintf(fo2,"%d %d\n", rows, cols);
    fprintf(fo2,"255\n");

    fprintf(fo3,"P5\n");
    fprintf(fo3,"%d %d\n", rows, cols);
    fprintf(fo3,"255\n");

    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    for (i=0;i<256;i++)
    { for (j=0;j<256;j++)
            {
              pic[i][j]  =  getc (fp1);
            }
    }

    for (p=-mr;p<=mr;p++)
    {  for (q=-mr;q<=mr;q++)
       {
          xmaskval = ((-q)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
          ymaskval = ((-p)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
          maskx[p+centy][q+centx] = xmaskval;
          masky[p+centy][q+centx] = ymaskval;
       }
    }

    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
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

    maxival = 0;
    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
      {
         mag[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                  (outpicy[i][j]*outpicy[i][j])));
         if (mag[i][j] > maxival)
            maxival = mag[i][j];

       }
    }

    for (i=0;i<256;i++)
      { for (j=0;j<256;j++)
        {
         mag[i][j] = (mag[i][j] / maxival) * 255;
         fprintf(fo1,"%c",(char)((int)(mag[i][j])));

        }
      }

    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
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

    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
         fprintf(fo2,"%c",(char)((int)(peak[i][j])));

        }
    }
    /*
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            if((peak[i][j] == 255.0) && (ival[i][j] > HI))
            {
                checkAdjacent(i, j);
            }
        }
    }

    for(i=0;i<256;i++)
    {
        printf("%f ", Histogram[i]);
    }*/

    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            Histogram[(int)mag[i][j]]++;
        }
    }

    CutOff = percent*rows*cols;

    for(i=255; i>=0;i--)
    {
        AreaOfTops += Histogram[i];
        if(AreaOfTops > CutOff)
            break;

    }

    HI=i;
    LO=(.35*HI);

    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            if(peak[i][j] == 255.0)
            {
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

    bool moretodo = true;
    while(moretodo == true)
    {
        moretodo =false;
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
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

    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
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


