//
//  main.c
//  TOV-final
//
//  Created by Weikang Lin on 5/30/18.
//  Copyright © 2018 Weikang Lin. All rights reserved.
//

# include<math.h>
# include<stdio.h>
# include<stdlib.h>
# define pi 3.1415
# define G 6.673e-11
# define c 2.99792458e8
# define M 1.99e30
# define dr 0.1

static double a[3000],b[3000];
static int j=0;

double che(double e);
double ch(double p);
double f(double r, double p, double e, double m);
double fm(double r,double e);

double che(double e)
{
    double le;
    double p,lp=0;
    int i=0;
    le=log10(e);
    i=j-1;
    for(;i>0;i=i-1)
    {
        if(le<=a[i]&&le>a[i-1])
        {
            lp=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
            break;
        }
    }
    p=pow(10,lp);
    return(p);
}

double ch(double p)
{
    double lp;
    double e,le=0;
    int i=0;
    lp=log10(p);
    i=j-1;
    for(;i>0;i=i-1)
    {
        if(lp<=b[i]&&lp>b[i-1])
        {
            le=a[i-1]*(lp-b[i])/(b[i-1]-b[i])+a[i]*(lp-b[i-1])/(b[i]-b[i-1]);
            break;
        }
    }
    e=pow(10,le);
    return(e);
}

double f(double r, double p, double e, double m)
{
    double g1, g2;
    g1=-G*(e+p/c/c)*(m+4*pi*r*r*r*p/c/c);
    g2=1/(r*r-2*G*r*m/c/c);
    return(g1*g2);
}

double fm(double r,double e)
{
    double mf;
    mf=4*pi*r*r*e;
    return(mf);
}



int main()
{
    FILE *inf,*outf;
    char filename[20],opfile[20];
    double r;
    double p, e, e0, m, y,m1,m2,m3,m4,p1,p2,p3,p4,power;
    printf("Please enter the data file's name: ");
    scanf("%s",filename);
    if((inf=fopen(filename,"r"))==NULL)
    {
        printf("file not found!");
        exit(0);
    }
    while(fscanf(inf,"%lf",&a[j])==1){fscanf(inf,"%lf",&b[j]);j++;}
    printf("Please enter a output file's name: ");
    scanf("%s",opfile);
    printf("processing...\n");
    outf=fopen(opfile,"w");
    power=pow(10.0,b[0]);
    e0=5E18;
    for(;e0>0.2E18;e0=pow(10.0,log10(e0)*0.999))
    {
        p=che(e0);
        m=4*dr*dr*pi*e0*dr/3;
        for(r=dr;p>power; r=r+dr)
        {
            e=ch(p);
            p1=f(r,p,e,m);
            m1=fm(r,e);
            if((p+dr*p1/2)>power) e=ch(p+dr*p1/2); else break;
            p2=f(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
            m2=fm(r+dr/2,e);
            if((p+dr*p2/2)>power) e=ch(p+dr*p2/2); else break;
            p3=f(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
            m3=fm(r+dr/2,e);
            if((p+dr*p3)>power) e=ch(p+dr*p3); else break;
            p4=f(r+dr,p+dr*p3,e,m+dr*m3);
            m4=fm(r+dr,e);
            p=p+dr*(p1+2*p2+2*p3+p4)/6;
            m=m+dr*(m1+2*m2+2*m3+m4)/6;
        }
        y=m/M;
        fprintf(outf,"%e",e0); /*-- Central density--*/
        fprintf(outf,"     ");
        fprintf(outf,"%f",y);  /*-- stellar mass --*/
        fprintf(outf,"     ");
        fprintf(outf,"%f",r);  /*-- stellar radius--*/
        fprintf(outf,"\n");
    }
    
    printf("Results have been saved in %s. \nPress any key to contiune\n",opfile);
    fclose(inf);
    fclose(outf);
    // getch();
}



