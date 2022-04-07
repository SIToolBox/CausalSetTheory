#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include<exception>
#include<iostream>


using namespace std;

#define N_CST_Store 43000
#define T_STEP 0  //if 1, timestep is 1Mpc time. if zero, it is 1e-4 in scale factor

// Random number generators
// ---------------------------------------
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

int lsprint = 0;
int error_index = 0;
//const int N_CST = 10000;
double N = 43000.0;
double alpha=0.01;  //controls the magnitude of the fluctuations
int maxi = 0; 


static double pi = 3.14159265;

double a_CST_Store[N_CST_Store + 1], tau_CST_Store[N_CST_Store + 1], time_CST_Store[N_CST_Store + 1], omegav_CST_Store[N_CST_Store + 1]; // Insert negative infinity
double vol_CST_Store[4][N_CST_Store + 1], omegaoH2_CST_Store[N_CST_Store + 1];
double omegav = 0.6925, toler = 1.0e-8;
int tstep = 0;  // if 1, time step is in t=1Mpc. if 0, updates every a=+1.0e-4
int wcount = 0; // counts number of times you hit negative H^2

// for saving the results to files
FILE *fHistory; //, *ftimeHistory, *ftauHistory, *fomegaH2History;


// According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.
// Returns a uniform random deviate between $0.0$ and $1.0$. Set idum to any negative value to initialize or reinitialize the sequence.
double ran3(long idum) 
{
    static int inext,inextp; 
    static long ma[56];        // The value 56 (range ma [1..55]) is special and
    static int iff=0;           // should not be modified; see Knuth.
    long mj, mk;
    int i, ii, k;
    if (idum <0 || iff == 0) { //Initialization.
        iff=1;
        mj=labs(MSEED-labs(idum));       //Initialize ma[55] using 
        mj %= MBIG;                       // large number MSEED. 
        ma[55]=mj;
        mk=1;
        for(i=1 ; i<=54 ; i++) {            // Now initialize the rest of the table,
            ii=(21 * i) % 55;               // in a slightly random order,
            ma[ii] = mk;                    // with numbers that are not especially random.
            mk = mj - mk;
            if(mk < MZ) mk += MBIG;
            mj = ma[ii];
            }
        for(k = 1; k <= 4; k++) // We randomize them by "warming up the generator.
            for ( i=1 ; i<=55 ; i++ ){
                ma[i] -= ma[1+(i+30) % 55];
                if(ma[i] < MZ) ma[i] +=MBIG;
                }
        //Prepare indices for our first generated number.
        //The constant 31 is special; see Knuth.
        inext=0;
        inextp=31;
        idum = 1;
        }
    //Here is where we start, except on initialization.
    if(++inext == 56) inext=1;          // Increment inext and inextp, wrapping around
    if (++inextp == 56) inextp =1 ;  // 56 to 1 .
    mj = ma[inext] - ma[inextp];  // Generate a new random number subtractively.
    if(mj < MZ) mj +=MBIG;
    ma[inext]=mj;
    return mj*FAC;
    }


class variables
{

public:
    double h0 = 70.74;

    double tcmb = 2.7254;
    double omegac = 0.2589;
    double omegab = 0.0486;
    double omegak = 0.0;
    double omegan = 0.0;
    int ndyn = 0;

    double grhom = (3.3379e-11) * h0 * h0;      // 3/c^2 = 3.3379e-11   ,this is 3H0^2/c^2 in Mpc^-2
    double grhog = (1.4952e-13) * pow(tcmb, 4); // 8*pi*G/(3c^2)*4*3*sigma_B/c^3 T^4           ,this is 8pi.G/(c^2) *(4sigma_B*T_present^4/c)/c^2 in Mpc^-2
    double grhor = (3.3957e-14) * pow(tcmb, 4); // 7/8*(4/11)^(4/3)*grhog (per neutrino species)
};

double rombint(double (*func)(double), double a, double b, double tol) // trapezoid approx to integral of func in interval [a,b]
{
    int MAXJ = 5;
    int MAXITER = 30, nint, i, j, k, jmax;
    double g[MAXJ + 1], h, gmax, g0, fourj, g1, error;

    h = 0.5 * (b - a);
    gmax = h * (func(a) + func(b)); // first trapezoid approx to integral of func in interval [a,b]
    g[0] = gmax;
    nint = 1;
    error = 1.e20;


    for (i = 0; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++)
    {
        /*     Calculate next trapezoidal rule approximation to integral. */
        g0 = 0.0;
        for (k = 1; k <= nint; k++)
            g0 += func(a + (k + k - 1) * h);

        g0 = 0.5 * g[0] + h * g0;
        h = 0.5 * h;
        nint *= 2;
        jmax = (i < MAXJ) ? i : MAXJ;
        fourj = 1.;

        for (j = 1; j <= jmax; j++)
        /* Use Richardson extrapolation. */
        {
            fourj *= 4.0;
            g1 = g0 + (g0 - g[j - 1]) / (fourj - 1.0);
            g[j - 1] = g0;
            g0 = g1;
        }

        /*  if (fabs(g0) > tol) */
        if (fabs(g0) > tol) // 1.e-30)
            error = 1.0 - gmax / g0;

        else
            error = gmax;
        gmax = g0;
        g[jmax] = g0;
    }

    if ((i > MAXITER) && (fabs(error) > tol))
        printf("rombint failed to converge; integral=%g, error=%g\n", g0, error);

    return g0;
}

double rombint2D(double (*func)(double, double), double c, double a, double b, double tol)
{
    
    int MAXJ = 5;
    int MAXITER = 30, nint, i, j, k, jmax;
    double g[MAXJ + 1], h, gmax, g0, fourj, g1, error;

    h = 0.5 * (b - a);
    gmax = h * (func(c, a) + func(c, b));
    g[0] = gmax;
    nint = 1;
    error = 1.e20;


    if(a<0.05) tol=1.0e-7;
    for (i = 0; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++)
    {
        /*     Calculate next trapezoidal rule approximation to integral. */
        g0 = 0.0;
        for (k = 1; k <= nint; k++)
            g0 += func(c, a + (k + k - 1) * h);

        g0 = 0.5 * g[0] + h * g0;
        h = 0.5 * h;
        nint *= 2;
        jmax = (i < MAXJ) ? i : MAXJ;
        fourj = 1.;

        for (j = 1; j <= jmax; j++)
        /* Use Richardson extrapolation. */
        {
            fourj *= 4.0;
            g1 = g0 + (g0 - g[j - 1]) / (fourj - 1.0);
            g[j - 1] = g0;
            g0 = g1;
        }

        /*  if (fabs(g0) > tol) */
        if (fabs(g0) > tol) // 1.e-30)
            error = 1.0 - gmax / g0;

        else
            error = gmax;
        gmax = g0;
        g[jmax] = g0;
    }

    if ((i > MAXITER) && (fabs(error) > tol))
        printf("rombint failed to converge; integral=%g, error=%g\n", g0, error);

    return g0;
}

double omegaDyn(double a)
{
    return omegav * pow(a, 4);
    }

double CSTdtauda(double a)
{
    variables Variables;
    double rx[2];
    double rhonu;
    double warned = 0;
    double dtauda1;
    double grho2; // Its a temporary variable for storing the total density

    grho2 = Variables.grhom * (Variables.omegac + Variables.omegab) * a + Variables.grhog + Variables.grhom * omegaDyn(a);

    if (grho2 < 0)
    {
        warned = 1; // hit a warning
        wcount = wcount + 1;
        throw grho2;
    }

    dtauda1 = sqrt(3.0 / grho2); 
    return (dtauda1);
}

double CSTdtimeda(double a)
{
    return a * CSTdtauda(a); 
}

double rand_gauss(double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double)X2);
    }

    do
    {
        U1 = -1 + ((double)ran3(10)) * 2;
        U2 = -1 + ((double)ran3(10)) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double)X1);
}

double CSTCalctau(double a)
{
    int i = (int)(a * N_CST_Store) + 1; // e-9 is just to ensure that you don't get into trouble with 2.9999...
    if(fabs(a*N_CST_Store - (int)(a * N_CST_Store)) < 1.0e-6)
        i = (int)(a * N_CST_Store-1.0e-6) + 1;    
    
    if(i>maxi) i=maxi;

    if (a > 1.0)
        i = N_CST_Store;
    if (a < 1.0e-8)
        return 0;
    if (i < 0)
    {
        tau_CST_Store[0] = 0.0;
        return 0.0;
    }
    if (tau_CST_Store[i] < 0.0)
    {
        if (tau_CST_Store[i - 1] < 0.0)
            CSTCalctau(a - 1.0 / N_CST_Store);
        tau_CST_Store[i] = tau_CST_Store[i - 1] + rombint(CSTdtauda, a_CST_Store[i - 1], a_CST_Store[i], toler);
    }

    if(lsprint==1)
        printf("\n(%e %e %e ) (%e %e %e %e %e[%d]) %e[%d] %e[%d]\n",a_CST_Store[i - 1], a, a_CST_Store[i], tau_CST_Store[i - 1], tau_CST_Store[i - 1] + rombint(CSTdtauda, a_CST_Store[i - 1], a, toler), tau_CST_Store[i - 1] + rombint(CSTdtauda, a_CST_Store[i - 1], a_CST_Store[i], toler), tau_CST_Store[i],tau_CST_Store[i+1],i,omegav_CST_Store[i-1],i-1,omegav_CST_Store[i],i);
    
    return tau_CST_Store[i - 1] + rombint(CSTdtauda, a_CST_Store[i - 1], a, toler);
}

double CSTCalcdVolumeda(double b, double a)
{
    constexpr static double pi = 3.14159265;
    double dV = (4 * pi / 3) * pow(a, 4) * pow((CSTCalctau(b) - CSTCalctau(a)), 3) * CSTdtauda(a); //
    return dV;
}

double CSTCalcdVolMomentda(double n, double a)
{
    return pow(a, 4) * pow(CSTCalctau(a), n) * CSTdtauda(a);
}

double CSTCalcVolume(double a)
{
    int j;
    int i = (int)(a * N_CST_Store) + 1;
    if(fabs(a*N_CST_Store - (int)(a * N_CST_Store)) < 1.0e-6)
        i = (int)(a * N_CST_Store - 1.0e-6) + 1;

    double mytest;
    if(i>maxi) i=maxi;
    
    if (a > 1.0)
        i = N_CST_Store;
    if (a < 1.0e-8)
        return 0;
    if (i <= 0)
    {
        for (j = 0; j < 4; j++)
            vol_CST_Store[j][0] = 0.0;
        return 0.0;
    }

    for (j = 0; j < 4; j++)
    {
        
        if (vol_CST_Store[j][i] < 0.0)
        {
            if (vol_CST_Store[j][i - 1] < 0)
                CSTCalcVolume(a - 1.0 / N_CST_Store);
            mytest = rombint2D(CSTCalcdVolMomentda, 1.0*j, a_CST_Store[i - 1], a_CST_Store[i], toler);
            vol_CST_Store[j][i] = vol_CST_Store[j][i - 1] + mytest;
            }
        }


    double tau = CSTCalctau(a), Vol = 0.0,Vol1 = 0.0,Vol2 = 0.0,Vol3 = 0.0;

    int factor[] = {1, -3, 3, -1};
    for (j = 0; j < 4; j++)
        Vol += factor[j] * pow(tau, 3 - j) * (vol_CST_Store[j][i - 1] + rombint2D(CSTCalcdVolMomentda, 1.0 * j, a_CST_Store[i - 1], a, toler));

    return (4 * pi / 3) * Vol;
}

void singlerun(int k)
{
    //k is the number of the run
    variables Variables;
    double a, taurmtest, taurmtest1;
    double aold, taurmtestold, taurmtest1old;
    double dVoldaold, Volold, dVolda, Vol;
    double dVol, S, Sold, rho;
    double zeta = 1.0;
    int j;
    double rhocr;
    double time1, timeold;
    double temp_time, temp_a, temp_aold;
    double presenttime = 1.0;
    double c2oG = 4.1553e+49; 
    omegav = 0.6925;

    //srand(time(0));
    aold = 1.0e-8;
    taurmtestold = 0.0;
    dVoldaold = 0.0;
    Volold = 0.0;
    S = 0.0;
    Sold = 0.0;
    rhocr = 2.7754e11 * 1.989e30 * pow(Variables.h0 / 100.0, 2); // critical density in kg/(Mpc)^3

    aold = 0.0;
    timeold = 0.0;
    temp_aold = 1.0e-8;

    for (int i = 0; i <= N_CST_Store; i++)
    {
        a_CST_Store[i] = i / N;
        omegav_CST_Store[i] = omegav;
        omegaoH2_CST_Store[i] = 0.0;
        time_CST_Store[i] = 0.0;
        tau_CST_Store[i] = -999999999.9; // Set everything initially to -infty
        for (j = 0; j < 4; j++)
            vol_CST_Store[j][i] = -999999999.9;
    }

    a_CST_Store[0] = 1.0e-8;
    tau_CST_Store[0] = 0.0;
    for (j = 0; j < 4; j++)
        vol_CST_Store[j][0] = 0.0;

    double H, omegaoH2;
    clock_t begin, end;

    for (int i = 0; i < N_CST_Store; i++)
    {
        a = i / N;
        omegav_CST_Store[i] = omegav;
        time1 = timeold + rombint(CSTdtimeda, aold, a, toler); // in steps of a, the codes integrates to find t(a)
        maxi = i;
        time_CST_Store[i] = time1;

        if (T_STEP == 1) // uses 1Mpc time steps 
        {

            if (time1 > presenttime && timeold < presenttime)
            {
                temp_a = aold + (presenttime - timeold) / (time1 - timeold) * (a - aold); //!!! is updating rho_lambda at temp_a,
                                                                                          // which has the same ratio of distance between a_old and a that presenttime has between timeold and time
                                                                                          // this is aimed to be the corresponding a(presenttime), but it's only approximate
                                                                                          // this is being used to integrate Friedmann between time and some newtime greater than presenttime+1;
                                                                                          // while we should use it to integrate from presenttime to presenttime+1
                presenttime += 1;

                Vol = CSTCalcVolume(temp_a);

                zeta = rand_gauss(0, 1);
                S = (Sold + alpha * zeta * sqrt(Vol - Volold));
                rho = S / Vol;
                
                //printf("%d %e %e %e %e %e %e %e\n", i, temp_a, presenttime, CSTCalctau(temp_a), CSTdtimeda(a), Vol, omegav, zeta);                
                fprintf(fHistory, "%d %e %e %e %e %e %e %e\n", i, temp_a, presenttime, CSTCalctau(temp_a), CSTdtimeda(a), Vol, omegav, zeta);
                fflush(stdout);

                // This rho is in (Mpc^-2)
                // We need to multiply it with c^2/G to convert it to kg . Mpc^-3
                // c = 299792458 m/s
                // G = 6.6743e-11 m^3 /kg /s^2
                // Mpc = 3.085678e22 m
                // c^2/G = ( c^2 / G ) 3.085678e22 = 4.1553e+49 kg/Mpc

                rho = rho * c2oG;
                omegav = rho / rhocr; // If its commented then it will give LCDM.

                Volold = Vol;
                Sold = S;
            }
        }
        else if (T_STEP == 0) // uses steps in a
        {
            if (i >= 30) // jump ahead the first 50 steps because the volume is too small, hence lambda is too big, hence the computation too slow
            {
                Vol = CSTCalcVolume(a);
                zeta = rand_gauss(0, 1);
                S = (Sold + alpha * zeta * sqrt(Vol - Volold));
                rho = S / Vol;

                // This rho is in (Mpc^-2)
                // We need to multiply it with c^2/G to convert it to kg . Mpc^-3
                // c = 299792458 m/s
                // G = 6.6743e-11 m^3 /kg /s^2
                // Mpc = 3.085678e22 m
                // c^2/G = ( c^2 / G ) 3.085678e22 = 4.1553e+49 kg/Mpc

                printf("%d %e %e %e %e %e\n", i, a, CSTdtimeda(a), Vol, omegav, zeta);
                fprintf(fHistory,"%d %e %e %e %e %e\n", i, a, CSTdtimeda(a), Vol, omegav, zeta);
                fflush(stdout);

                rho = rho * c2oG;

                omegav = rho / rhocr; // If its commented then it will give LCDM.
                Volold = Vol;
                Sold = S;
            }
        }
        else if (T_STEP == 2) // uses volume of past light cone as time coordinate
        {
            Vol = CSTCalcVolume(a);
            if (Vol - Volold > 1) // only update omegav if .1 (Mpc^4) of volume has been increased
            {

                zeta = rand_gauss(0, 1);
                S = (Sold + alpha * zeta * sqrt(Vol - Volold));
                if (Vol <= 0.0) //#to avoid infinity at first step
                {
                    rho = 0;
                }
                else
                {
                    rho = S / Vol;
                }
                rho = rho * c2oG;
                omegav = rho / rhocr; // If its commented then it will give LCDM.
                Volold = Vol;
                Sold = S;
            }
        }

        aold = a;
        timeold = time1;
        while(presenttime<timeold)  
            presenttime++;
    }

    printf("\n Step: %d, Terminations:%d \n", k, wcount);
}

double calcChi2(double omegaH2, double alpha, int k)
{
    try{
        
    }
    catch(...)
    {
        return 999999999.0;
    }
}

int mcmc()
{
    // Initial 
    double omegamH2, alpha, old_omegamH2, old_alpha, delta_omegamH2_alpha, delta_alpha_omegamH2, delta_omegamH2, delta_alpha; 
    double low_omegaH2 = 0.12, high_omegaH2 = 0.16, std_omegaH2 = 0.003;
    double low_alpha = 0.005, high_alpha = 0.025, std_alpha = 0.003;

    

    for( int jj = 0; jj<50000; jj++)
    {
        do
        {
            delta_omegamH2_alpha = rand_gauss(0.0,1.0); 
            delta_omegamH2 = delta_omegamH2_alpha *  std_omegaH2;
            omegamH2 = old_omegamH2 + delta_omegamH2;
            } while(( omegaH2 < low_omegaH2 ) || ( omegaH2 > high_omegaH2 ));

        do
        {
            delta_alpha_omegamH2 = rand_gauss(0.0,1.0);
            delta_alpha = delta_alpha_omegamH2 * std_alpha;
            alpha = old_alpha + delta_alpha;
            } while(( alpha < low_alpha ) || ( alpha > high_alpha ));
  
        Chi2 = calcChi2(omegaH2, alpha, jj);

        if( Chi2 < Chi2old )    // Accept It
        {
            old_alpha = alpha;
            old_omegamH2 = omegaH2;
            Chi2old = Chi2;
            }
        else
        {
            u = exp(Chi2old - Chi2);
            if(u > ran3())
            {
                old_alpha = alpha;
                old_omegamH2 = omegaH2;
                Chi2old = Chi2;
                }
            }
        }

    printf("%d %e %e %e \n", jj, old_omegamH2, old_alpha, 2*Chi2);
    }



int main()
{
    char sHistory[500], command[500];
    int srandom;
    int startstep, stepnumber;
    printf("\nEnter Start step and Iteration per step:\n");
    scanf("%d %d",&startstep, &stepnumber);

    for (int i = startstep*stepnumber; i < (startstep+1)*stepnumber; i++) //R
    {
        sprintf(sHistory,"History_%d.d",i);
        fHistory = fopen(sHistory,"w");
        try{
            srandom = -1*i;
            ran3(srandom);
            singlerun(i);
            fclose(fHistory);
        }
        catch(...)
        {
            printf("\n Negative H^2 occured\n");
            fclose(fHistory);
            sprintf(command,"N_History_%d.d",i);
            rename(sHistory, command);
        }

    } // end of the for loop for repeating the run

    printf("\n Number of terminations:%d \n", wcount);
}
