#include "math.h"
#include<stdio.h>

/*  Merging the nlplant.c (lofi) and nlplant_hifi.c to use
    same equations of motion, navigation equations and use 
    own look-up tables decided by a flag.                   */
    
void atmos(double,double,double*);          /* Used by both */
void accels(double*,double*,double*);       /* Used by both */

#include "lofi_F16_AeroData.c"              /* LOFI Look-up header file*/
#include "hifi_F16_AeroData.c"              /* HIFI Look-up header file*/

void Nlplant(double*,double*,int);

/*########################################*/
/*### Added for mex function in matlab ###*/
/*########################################*/

/*########################################*/
/*########################################*/

void Nlplant(double *xu, double *xdot, int fidelity){ //I modified to have fidelity as a seperate input!

  int fi_flag;

  /* #include f16_constants */
  double g    = 32.17;          /* gravity, ft/s^2 */
  double m    = 636.94;         /* mass, slugs */
  double B    = 30.0;             /* span, ft */
  double S    = 300.0;            /* planform area, ft^2 */
  double cbar = 11.32;          /* mean aero chord, ft */
  double xcgr = 0.35;      /* reference center of gravity as a fraction of cbar */
  double xcg  = 0.25;      /* center of gravity as a fraction of cbar. */

  double Heng = 0.0;              /* turbine momentum along roll axis. */
  double pi   = acos(-1);
  double r2d;                   /* radians to degrees */


/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double Jy  = 55814.0;           /* slug-ft^2 */ 
  double Jxz = 982.0;             /* slug-ft^2 */     
  double Jz  = 63100.0;           /* slug-ft^2 */
  double Jx  = 9496.0;            /* slug-ft^2 */

  double *temp;

  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double T, el, ail, rud, dail, drud, lef, dlef;
  double qbar, mach, ps;
  double U, V, W, Udot,Vdot,Wdot;
  double L_tot, M_tot, N_tot, denom;

  double Cx_tot, Cx, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm, eta_el, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_ds;
  double Cy_tot, Cy, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;

  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/

  r2d  = 180.0/pi;     /* radians to degrees */

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */


  npos  = xu[0];   /* north position */
  epos  = xu[1];   /* east position */
  alt   = xu[2];   /* altitude */
  phi   = xu[3];   /* orientation angles in rad. */
  theta = xu[4];
  psi   = xu[5];

  vt    = xu[6];     /* total velocity */
  alpha = xu[7]*r2d; /* angle of attack in degrees */
  beta  = xu[8]*r2d; /* sideslip angle in degrees */
  P     = xu[9];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = xu[10];    /* Pitch Rate--- pitching moment is M */
  R     = xu[11];    /* Yaw Rate  --- yawing   moment is N */

  sa    = sin(xu[7]); /* sin(alpha) */
  ca    = cos(xu[7]); /* cos(alpha) */
  sb    = sin(xu[8]); /* sin(beta)  */
  cb    = cos(xu[8]); /* cos(beta)  */
  tb    = tan(xu[8]); /* tan(beta)  */

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  T     = xu[12];   /* thrust */
  el    = xu[13];   /* Elevator setting in degrees. */
  ail   = xu[14];   /* Ailerons mex setting in degrees. */
  rud   = xu[15];   /* Rudder setting in degrees. */
  lef   = xu[16];   /* Leading edge flap setting in degrees */
  
  //fi_flag = xu[17]/1;       /* fi_flag */
  fi_flag = fidelity; // Johns change to pass in integer seperate from numpy float array
    
  /* dail  = ail/20.0;   aileron normalized against max angle */
  /* The aileron was normalized using 20.0 but the NASA report and
     S&L both have 21.5 deg. as maximum deflection. */
  /* As a result... */
  dail  = ail/21.5;
  drud  = rud/30.0;  /* rudder normalized against max angle */
  dlef  = (1 - lef/25.0);  /* leading edge flap normalized against max angle */


  /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure and mach number
     %%%%%%%%%%%%%%%%%% */

atmos(alt,vt,temp);
   mach = temp[0];
   qbar = temp[1];
   ps   = temp[2];

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  /* %%%%%%%%%%%%%%%%%%
     Navigation Equations
     %%%%%%%%%%%%%%%%%% */

   U = vt*ca*cb;  /* directional velocities. */
   V = vt*sb;
   W = vt*sa*cb;

/* nposdot */
xdot[0] = U*(ct*cpsi) + 
            V*(sphi*cpsi*st - cphi*spsi) + 
            W*(cphi*st*cpsi + sphi*spsi);

/* eposdot */  
xdot[1] = U*(ct*spsi) + 
            V*(sphi*spsi*st + cphi*cpsi) + 
            W*(cphi*st*spsi - sphi*cpsi);

/* altdot */
xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct);

  /* %%%%%%%%%%%%%%%%%%%
     Kinematic equations
     %%%%%%%%%%%%%%%%%%% */
/* phidot */
xdot[3] = P + tt*(Q*sphi + R*cphi);


/* theta dot */
xdot[4] = Q*cphi - R*sphi;

/* psidot */
xdot[5] = (Q*sphi + R*cphi)/ct;


/* %%%%%%%%%%%%%%%%%%
        Table lookup
     %%%%%%%%%%%%%%%%%% */

if (fi_flag == 1)          /* HIFI Table */
{
    hifi_C(alpha,beta,el,temp);
        Cx = temp[0];
        Cz = temp[1];
        Cm = temp[2];
        Cy = temp[3];
        Cn = temp[4];
        Cl = temp[5];

    hifi_damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    hifi_C_lef(alpha,beta,temp);
        delta_Cx_lef = temp[0];
        delta_Cz_lef = temp[1];
        delta_Cm_lef = temp[2];
        delta_Cy_lef = temp[3];
        delta_Cn_lef = temp[4];
        delta_Cl_lef = temp[5];

    hifi_damping_lef(alpha,temp);
        delta_Cxq_lef = temp[0];
        delta_Cyr_lef = temp[1];
        delta_Cyp_lef = temp[2];
        delta_Czq_lef = temp[3];
        delta_Clr_lef = temp[4];
        delta_Clp_lef = temp[5];
        delta_Cmq_lef = temp[6];
        delta_Cnr_lef = temp[7];
        delta_Cnp_lef = temp[8];

    hifi_rudder(alpha,beta,temp);
        delta_Cy_r30 = temp[0];
        delta_Cn_r30 = temp[1];
        delta_Cl_r30 = temp[2];

    hifi_ailerons(alpha,beta,temp);
        delta_Cy_a20     = temp[0];
        delta_Cy_a20_lef = temp[1];
        delta_Cn_a20     = temp[2];
        delta_Cn_a20_lef = temp[3];
        delta_Cl_a20     = temp[4];
        delta_Cl_a20_lef = temp[5];

    hifi_other_coeffs(alpha,el,temp);
        delta_Cnbeta = temp[0];
        delta_Clbeta = temp[1];
        delta_Cm     = temp[2];
        eta_el       = temp[3];
        delta_Cm_ds  = 0;        /* ignore deep-stall effect */
    
}

else if (fi_flag == 0)
{     
/* ##############################################
   ##########LOFI Table Look-up #################
   ##############################################*/

/* The lofi model does not include the
   leading edge flap.  All terms multiplied
   dlef have been set to zero but just to 
   be sure we will set it to zero. */
    
    dlef = 0.0;     

    damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    dmomdcon(alpha,beta, temp);
        delta_Cl_a20 = temp[0];     /* Formerly dLda in nlplant.c */
        delta_Cl_r30 = temp[1];     /* Formerly dLdr in nlplant.c */
        delta_Cn_a20 = temp[2];     /* Formerly dNda in nlplant.c */
        delta_Cn_r30 = temp[3];     /* Formerly dNdr in nlplant.c */

    clcn(alpha,beta,temp);
        Cl = temp[0];
        Cn = temp[1];

    cxcm(alpha,el,temp);
        Cx = temp[0];
        Cm = temp[1];

    Cy = -.02*beta + .021*dail + .086*drud;

    cz(alpha,beta,el,temp);
        Cz = temp[0];
/*##################################################
        
        
/*##################################################
  ##  Set all higher order terms of hifi that are ##
  ##  not applicable to lofi equal to zero. ########
  ##################################################*/
     
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0;     /* Needs to be one. See equation for Cm_tot*/
        delta_Cm_ds     = 0.0;
                     
/*##################################################
  ##################################################*/ 
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
(as on NASA report p37-40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* XXXXXXXX Cx_tot XXXXXXXX */

dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);

Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

   /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 

dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);

Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

   /* MMMMMMMM Cm_tot MMMMMMMM */ 

dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);

Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

   /* YYYYYYYY Cy_tot YYYYYYYY */

dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;

dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);

dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);

Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

   /* NNNNNNNN Cn_tot NNNNNNNN */ 

dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;

dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);

dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);

Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

   /* LLLLLLLL Cl_tot LLLLLLLL */

dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;

dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);

dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);

Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      compute Udot,Vdot, Wdot,(as on NASA report p36)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vt_dot equation (from S&L, p82)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;

   /* %%%%%%%%%%%%%%%%%%
      alpha_dot equation
      %%%%%%%%%%%%%%%%%% */

xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);

  /* %%%%%%%%%%%%%%%%%
     beta_dot equation
     %%%%%%%%%%%%%%%%% */

xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);



  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Pdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[9] =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;


  /* %%%%%%%%%%%%%%%%%%%%%%%
     Qdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Rdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

/*########################################*/
/*### Create accelerations anx_cg, any_cg */
/*### ans anz_cg as outputs ##############*/
/*########################################*/

accels(xu,xdot,temp);

xdot[12]  = temp[0];	/* anx_cg */
xdot[13]  = temp[1];	/* any_cg */
xdot[14]  = temp[2];	/* anz_cg */
xdot[15]  = mach;
xdot[16]  = qbar;
xdot[17]  = ps;

/*########################################*/
/*########################################*/

free(temp);

}; /*##### END of nlplant() ####*/

/*########################################*/
/*### Called Sub-Functions  ##############*/
/*########################################*/

/*########################################*/
/* Function for mach and qbar */
/*########################################*/

void atmos(double alt, double vt, double *coeff ){

    double rho0 = 2.377e-3;
    double tfac, temp, rho, mach, qbar, ps;

    tfac =1 - .703e-5*(alt);
    temp = 519.0*tfac;
    if (alt >= 35000.0) {
       temp=390;
    }

    rho=rho0*pow(tfac,4.14);
    mach = (vt)/sqrt(1.4*1716.3*temp);
    qbar = .5*rho*pow(vt,2);
    ps   = 1715.0*rho*temp;

    if (ps == 0){
      ps = 1715;
      }

    coeff[0] = mach;
    coeff[1] = qbar;
    coeff[2] = ps;
}

/*########################################*/
/*########################################*/


/*########################################*/
/*### Port from matlab fix() function ####*/
/*########################################*/

/* port from matlab sign() function */


/*########################################*/
/*########################################*/


/*########################################*/
/*### Calculate accelerations from states */
/*### and state derivatives. ############ */
/*########################################*/

void accels(double *state,
            double *xdot,
            double *y)

{

#define grav 32.174 

double sina, cosa, sinb, cosb ;
double vel_u, vel_v, vel_w ;
double u_dot, v_dot, w_dot ;
double nx_cg, ny_cg, nz_cg ;

sina = sin(state[7]) ;
cosa = cos(state[7]) ;
sinb = sin(state[8]) ;
cosb = cos(state[8]) ;
vel_u = state[6]*cosb*cosa ;
vel_v = state[6]*sinb ;
vel_w = state[6]*cosb*sina ;
u_dot =          cosb*cosa*xdot[6]
      - state[6]*sinb*cosa*xdot[8] 
      - state[6]*cosb*sina*xdot[7] ;
v_dot =          sinb*xdot[6] 
      + state[6]*cosb*xdot[8] ;
w_dot =          cosb*sina*xdot[6]
      - state[6]*sinb*sina*xdot[8] 
      + state[6]*cosb*cosa*xdot[7] ;
nx_cg = 1.0/grav*(u_dot + state[10]*vel_w - state[11]*vel_v)
      + sin(state[4]) ;
ny_cg = 1.0/grav*(v_dot + state[11]*vel_u - state[9]*vel_w)
      - cos(state[4])*sin(state[3]) ;
nz_cg = -1.0/grav*(w_dot + state[9]*vel_v - state[10]*vel_u)
      + cos(state[4])*cos(state[3]) ;

y[0] = nx_cg ;
y[1] = ny_cg ;
y[2] = nz_cg ;


} 

/*########################################*/
/*########################################*/

/*########################################*/
/*########################################*/

void Jac(double *xu, double *xdot, int fidelity, double *jac){ //I modified to have fidelity as a seperate input!

  int fi_flag;

  /* #include f16_constants */
  double g    = 32.17;          /* gravity, ft/s^2 */
  double m    = 636.94;         /* mass, slugs */
  double B    = 30.0;             /* span, ft */
  double S    = 300.0;            /* planform area, ft^2 */
  double cbar = 11.32;          /* mean aero chord, ft */
  double xcgr = 0.25;      /* reference center of gravity as a fraction of cbar */
  double xcg  = 0.25;      /* center of gravity as a fraction of cbar. */

  double Heng = 0.0;              /* turbine momentum along roll axis. */
  double pi   = acos(-1);
  double r2d;                   /* radians to degrees */


/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double Jy  = 55814.0;           /* slug-ft^2 */ 
  double Jxz = 982.0;             /* slug-ft^2 */     
  double Jz  = 63100.0;           /* slug-ft^2 */
  double Jx  = 9496.0;            /* slug-ft^2 */

  double *temp;

  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double T, el, ail, rud, dail, drud, lef, dlef;
  double qbar, mach, ps;
  double U, V, W, Udot,Vdot,Wdot;
  double L_tot, M_tot, N_tot, denom;

  double Cx_tot, Cx, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm, eta_el, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_ds;
  double Cy_tot, Cy, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;

  double dx0dot_dx0, dx0dot_dx1, dx0dot_dx2, dx0dot_dx3, dx0dot_dx4, dx0dot_dx5, dx0dot_dx6, dx0dot_dx7, dx0dot_dx8, dx0dot_dx9, dx0dot_dx10, dx0dot_dx11;
  double dx1dot_dx0, dx1dot_dx1, dx1dot_dx2, dx1dot_dx3, dx1dot_dx4, dx1dot_dx5, dx1dot_dx6, dx1dot_dx7, dx1dot_dx8, dx1dot_dx9, dx1dot_dx10, dx1dot_dx11;
  double dx2dot_dx0, dx2dot_dx1, dx2dot_dx2, dx2dot_dx3, dx2dot_dx4, dx2dot_dx5, dx2dot_dx6, dx2dot_dx7, dx2dot_dx8, dx2dot_dx9, dx2dot_dx10, dx2dot_dx11;
  double dx3dot_dx0, dx3dot_dx1, dx3dot_dx2, dx3dot_dx3, dx3dot_dx4, dx3dot_dx5, dx3dot_dx6, dx3dot_dx7, dx3dot_dx8, dx3dot_dx9, dx3dot_dx10, dx3dot_dx11;
  double dx4dot_dx0, dx4dot_dx1, dx4dot_dx2, dx4dot_dx3, dx4dot_dx4, dx4dot_dx5, dx4dot_dx6, dx4dot_dx7, dx4dot_dx8, dx4dot_dx9, dx4dot_dx10, dx4dot_dx11;
  double dx5dot_dx0, dx5dot_dx1, dx5dot_dx2, dx5dot_dx3, dx5dot_dx4, dx5dot_dx5, dx5dot_dx6, dx5dot_dx7, dx5dot_dx8, dx5dot_dx9, dx5dot_dx10, dx5dot_dx11;
  double dx6dot_dx0, dx6dot_dx1, dx6dot_dx2, dx6dot_dx3, dx6dot_dx4, dx6dot_dx5, dx6dot_dx6, dx6dot_dx7, dx6dot_dx8, dx6dot_dx9, dx6dot_dx10, dx6dot_dx11;
  double dx7dot_dx0, dx7dot_dx1, dx7dot_dx2, dx7dot_dx3, dx7dot_dx4, dx7dot_dx5, dx7dot_dx6, dx7dot_dx7, dx7dot_dx8, dx7dot_dx9, dx7dot_dx10, dx7dot_dx11;
  double dx8dot_dx0, dx8dot_dx1, dx8dot_dx2, dx8dot_dx3, dx8dot_dx4, dx8dot_dx5, dx8dot_dx6, dx8dot_dx7, dx8dot_dx8, dx8dot_dx9, dx8dot_dx10, dx8dot_dx11;
  double dx9dot_dx0, dx9dot_dx1, dx9dot_dx2, dx9dot_dx3, dx9dot_dx4, dx9dot_dx5, dx9dot_dx6, dx9dot_dx7, dx9dot_dx8, dx9dot_dx9, dx9dot_dx10, dx9dot_dx11;
  double dx10dot_dx0, dx10dot_dx1, dx10dot_dx2, dx10dot_dx3, dx10dot_dx4, dx10dot_dx5, dx10dot_dx6, dx10dot_dx7, dx10dot_dx8, dx10dot_dx9, dx10dot_dx10, dx10dot_dx11;
  double dx11dot_dx0, dx11dot_dx1, dx11dot_dx2, dx11dot_dx3, dx11dot_dx4, dx11dot_dx5, dx11dot_dx6, dx11dot_dx7, dx11dot_dx8, dx11dot_dx9, dx11dot_dx10, dx11dot_dx11;

  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/

  r2d  = 180.0/pi;     /* radians to degrees */

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */


  npos  = xu[0];   /* north position */
  epos  = xu[1];   /* east position */
  alt   = xu[2];   /* altitude */
  phi   = xu[3];   /* orientation angles in rad. */
  theta = xu[4];
  psi   = xu[5];

  vt    = xu[6];     /* total velocity */
  alpha = xu[7]*r2d; /* angle of attack in degrees */
  beta  = xu[8]*r2d; /* sideslip angle in degrees */
  P     = xu[9];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = xu[10];    /* Pitch Rate--- pitching moment is M */
  R     = xu[11];    /* Yaw Rate  --- yawing   moment is N */

  sa    = sin(xu[7]); /* sin(alpha) */
  ca    = cos(xu[7]); /* cos(alpha) */
  sb    = sin(xu[8]); /* sin(beta)  */
  cb    = cos(xu[8]); /* cos(beta)  */
  tb    = tan(xu[8]); /* tan(beta)  */

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  T     = xu[12];   /* thrust */
  el    = xu[13];   /* Elevator setting in degrees. */
  ail   = xu[14];   /* Ailerons mex setting in degrees. */
  rud   = xu[15];   /* Rudder setting in degrees. */
  lef   = xu[16];   /* Leading edge flap setting in degrees */
  
  //fi_flag = xu[17]/1;       /* fi_flag */
  fi_flag = fidelity; // Johns change to pass in integer seperate from numpy float array
    
  /* dail  = ail/20.0;   aileron normalized against max angle */
  /* The aileron was normalized using 20.0 but the NASA report and
     S&L both have 21.5 deg. as maximum deflection. */
  /* As a result... */
  dail  = ail/21.5;
  drud  = rud/30.0;  /* rudder normalized against max angle */
  dlef  = (1 - lef/25.0);  /* leading edge flap normalized against max angle */


  /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure and mach number
     %%%%%%%%%%%%%%%%%% */

atmos(alt,vt,temp);
   mach = temp[0];
   qbar = temp[1];
   ps   = temp[2];

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  /* %%%%%%%%%%%%%%%%%%
     Navigation Equations
     %%%%%%%%%%%%%%%%%% */

   U = vt*ca*cb;  /* directional velocities. */
   V = vt*sb;
   W = vt*sa*cb;

/* nposdot */
xdot[0] = U*(ct*cpsi) + 
            V*(sphi*cpsi*st - cphi*spsi) + 
            W*(cphi*st*cpsi + sphi*spsi);

/* eposdot */  
xdot[1] = U*(ct*spsi) + 
            V*(sphi*spsi*st + cphi*cpsi) + 
            W*(cphi*st*spsi - sphi*cpsi);

/* altdot */
xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct);

  /* %%%%%%%%%%%%%%%%%%%
     Kinematic equations
     %%%%%%%%%%%%%%%%%%% */
/* phidot */
xdot[3] = P + tt*(Q*sphi + R*cphi);


/* theta dot */
xdot[4] = Q*cphi - R*sphi;

/* psidot */
xdot[5] = (Q*sphi + R*cphi)/ct;


/* %%%%%%%%%%%%%%%%%%
        Table lookup
     %%%%%%%%%%%%%%%%%% */

if (fi_flag == 1)          /* HIFI Table */
{
    hifi_C(alpha,beta,el,temp);
        Cx = temp[0];
        Cz = temp[1];
        Cm = temp[2];
        Cy = temp[3];
        Cn = temp[4];
        Cl = temp[5];

    hifi_damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    hifi_C_lef(alpha,beta,temp);
        delta_Cx_lef = temp[0];
        delta_Cz_lef = temp[1];
        delta_Cm_lef = temp[2];
        delta_Cy_lef = temp[3];
        delta_Cn_lef = temp[4];
        delta_Cl_lef = temp[5];

    hifi_damping_lef(alpha,temp);
        delta_Cxq_lef = temp[0];
        delta_Cyr_lef = temp[1];
        delta_Cyp_lef = temp[2];
        delta_Czq_lef = temp[3];
        delta_Clr_lef = temp[4];
        delta_Clp_lef = temp[5];
        delta_Cmq_lef = temp[6];
        delta_Cnr_lef = temp[7];
        delta_Cnp_lef = temp[8];

    hifi_rudder(alpha,beta,temp);
        delta_Cy_r30 = temp[0];
        delta_Cn_r30 = temp[1];
        delta_Cl_r30 = temp[2];

    hifi_ailerons(alpha,beta,temp);
        delta_Cy_a20     = temp[0];
        delta_Cy_a20_lef = temp[1];
        delta_Cn_a20     = temp[2];
        delta_Cn_a20_lef = temp[3];
        delta_Cl_a20     = temp[4];
        delta_Cl_a20_lef = temp[5];

    hifi_other_coeffs(alpha,el,temp);
        delta_Cnbeta = temp[0];
        delta_Clbeta = temp[1];
        delta_Cm     = temp[2];
        eta_el       = temp[3];
        delta_Cm_ds  = 0;        /* ignore deep-stall effect */
    
}

else if (fi_flag == 0)
{     
/* ##############################################
   ##########LOFI Table Look-up #################
   ##############################################*/

/* The lofi model does not include the
   leading edge flap.  All terms multiplied
   dlef have been set to zero but just to 
   be sure we will set it to zero. */
    
    dlef = 0.0;     

    damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    dmomdcon(alpha,beta, temp);
        delta_Cl_a20 = temp[0];     /* Formerly dLda in nlplant.c */
        delta_Cl_r30 = temp[1];     /* Formerly dLdr in nlplant.c */
        delta_Cn_a20 = temp[2];     /* Formerly dNda in nlplant.c */
        delta_Cn_r30 = temp[3];     /* Formerly dNdr in nlplant.c */

    clcn(alpha,beta,temp);
        Cl = temp[0];
        Cn = temp[1];

    cxcm(alpha,el,temp);
        Cx = temp[0];
        Cm = temp[1];

    Cy = -.02*beta + .021*dail + .086*drud;

    cz(alpha,beta,el,temp);
        Cz = temp[0];
/*##################################################
        
        
/*##################################################
  ##  Set all higher order terms of hifi that are ##
  ##  not applicable to lofi equal to zero. ########
  ##################################################*/
     
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0;     /* Needs to be one. See equation for Cm_tot*/
        delta_Cm_ds     = 0.0;
                     
/*##################################################
  ##################################################*/ 
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
(as on NASA report p37-40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* XXXXXXXX Cx_tot XXXXXXXX */

dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);

Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

   /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 

dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);

Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

   /* MMMMMMMM Cm_tot MMMMMMMM */ 

dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);

Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

   /* YYYYYYYY Cy_tot YYYYYYYY */

dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;

dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);

dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);

Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

   /* NNNNNNNN Cn_tot NNNNNNNN */ 

dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;

dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);

dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);

Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

   /* LLLLLLLL Cl_tot LLLLLLLL */

dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;

dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);

dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);

Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      compute Udot,Vdot, Wdot,(as on NASA report p36)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vt_dot equation (from S&L, p82)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;

   /* %%%%%%%%%%%%%%%%%%
      alpha_dot equation
      %%%%%%%%%%%%%%%%%% */

xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);

  /* %%%%%%%%%%%%%%%%%
     beta_dot equation
     %%%%%%%%%%%%%%%%% */

xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);



  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Pdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[9] =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;


  /* %%%%%%%%%%%%%%%%%%%%%%%
     Qdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Rdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

/*########################################*/
/*### Create accelerations anx_cg, any_cg */
/*### ans anz_cg as outputs ##############*/
/*########################################*/

accels(xu,xdot,temp);

xdot[12]  = temp[0];	/* anx_cg */
xdot[13]  = temp[1];	/* any_cg */
xdot[14]  = temp[2];	/* anz_cg */
xdot[15]  = mach;
xdot[16]  = qbar;
xdot[17]  = ps;

/*########################################*/
/*########################################*/

// npos_dot s relationship with other primary states
dx0dot_dx0 = 0;
dx0dot_dx1 = 0;
dx0dot_dx2 = 0; 
dx0dot_dx3 = sb*vt*(sphi*spsi + cphi*cpsi*st) + cb*sa*vt*(cphi*spsi - cpsi*sphi*st);
dx0dot_dx4 = cpsi*ct*sb*sphi*vt - ca*cb*cpsi*st*vt + cb*cphi*cpsi*ct*sa*vt;
dx0dot_dx5 = cb*sa*vt*(cpsi*sphi - cphi*spsi*st) - sb*vt*(cphi*cpsi + sphi*spsi*st) - ca*cb*ct*spsi*vt;
dx0dot_dx6 = cb*sa*(sphi*spsi + cphi*cpsi*st) - sb*(cphi*spsi - cpsi*sphi*st) + ca*cb*cpsi*ct;
dx0dot_dx7 = ca*cb*vt*(sphi*spsi + cphi*cpsi*st) - cb*cpsi*ct*sa*vt;
dx0dot_dx8 = - cb*vt*(cphi*spsi - cpsi*sphi*st) - sa*sb*vt*(sphi*spsi + cphi*cpsi*st) - ca*cpsi*ct*sb*vt;
dx0dot_dx9 = 0; 
dx0dot_dx10 = 0; 
dx0dot_dx11 = 0;

// epos_dot s relationship with other primary states
dx1dot_dx0 = 0;
dx1dot_dx1 = 0;
dx1dot_dx2 = 0;
dx1dot_dx3 = - sb*vt*(cpsi*sphi - cphi*spsi*st) - cb*sa*vt*(cphi*cpsi + sphi*spsi*st);
dx1dot_dx4 = ct*sb*sphi*spsi*vt - ca*cb*spsi*st*vt + cb*cphi*ct*sa*spsi*vt;
dx1dot_dx5 = cb*sa*vt*(sphi*spsi + cphi*cpsi*st) - sb*vt*(cphi*spsi - cpsi*sphi*st) + ca*cb*cpsi*ct*vt;
dx1dot_dx6 = sb*(cphi*cpsi + sphi*spsi*st) - cb*sa*(cpsi*sphi - cphi*spsi*st) + ca*cb*ct*spsi;
dx1dot_dx7 = - ca*cb*vt*(cpsi*sphi - cphi*spsi*st) - cb*ct*sa*spsi*vt;
dx1dot_dx8 = cb*vt*(cphi*cpsi + sphi*spsi*st) + sa*sb*vt*(cpsi*sphi - cphi*spsi*st) - ca*ct*sb*spsi*vt;
dx1dot_dx9 = 0;
dx1dot_dx10 = 0;
dx1dot_dx11 = 0;

// h_dot s relationship with other primary states
dx2dot_dx0 = 0;
dx2dot_dx1 = 0;
dx2dot_dx2 = 0;
dx2dot_dx3 = cb*ct*sa*sphi*vt - cphi*ct*sb*vt;
dx2dot_dx4 = ca*cb*ct*vt + sb*sphi*st*vt + cb*cphi*sa*st*vt;
dx2dot_dx5 = 0;
dx2dot_dx6 = ca*cb*st - ct*sb*sphi - cb*cphi*ct*sa;
dx2dot_dx7 = - cb*sa*st*vt - ca*cb*cphi*ct*vt;
dx2dot_dx8 = cphi*ct*sa*sb*vt - ca*sb*st*vt - cb*ct*sphi*vt;
dx2dot_dx9 = 0;
dx2dot_dx10 = 0;
dx2dot_dx11 = 0;

// phi_dot s relationship with other primary states
dx3dot_dx0 = 0;
dx3dot_dx1 = 0;
dx3dot_dx2 = 0;
dx3dot_dx3 = tt*(Q*cphi - R*sphi);
dx3dot_dx4 = (tt*tt + 1)*(R*cphi + Q*sphi);
dx3dot_dx5 = 0;
dx3dot_dx6 = 0;
dx3dot_dx7 = 0;
dx3dot_dx8 = 0;
dx3dot_dx9 = 1;
dx3dot_dx10 = sphi*tt;
dx3dot_dx11 = cphi*tt;

// theta_dot s relationship with other primary states
dx4dot_dx0 = 0;
dx4dot_dx1 = 0;
dx4dot_dx2 = 0;
dx4dot_dx3 = - R*cphi - Q*sphi;
dx4dot_dx4 = 0;
dx4dot_dx5 = 0;
dx4dot_dx6 = 0;
dx4dot_dx7 = 0;
dx4dot_dx8 = 0;
dx4dot_dx9 = 0;
dx4dot_dx10 = cphi;
dx4dot_dx11 = -sphi;

// psi_dot s relationship with other primary states
dx5dot_dx0 = 0;
dx5dot_dx1 = 0;
dx5dot_dx2 = 0;
dx5dot_dx3 = (Q*cphi - R*sphi)/ct;
dx5dot_dx4 = (st*(R*cphi + Q*sphi))/(ct*ct);
dx5dot_dx5 = 0;
dx5dot_dx6 = 0;
dx5dot_dx7 = 0;
dx5dot_dx8 = 0;
dx5dot_dx9 = 0;
dx5dot_dx10 = sphi/ct;
dx5dot_dx11 = cphi/ct;

// vt_dot s relationship with other primary states
dx6dot_dx0 = 0;
dx6dot_dx1 = 0;
dx6dot_dx2 = 0;
dx6dot_dx3 = ((3217*cphi*ct*sb*vt)/100 - (3217*cb*ct*sa*sphi*vt)/100)/vt;
dx6dot_dx4 = -((3217*ca*cb*ct*vt)/100 + (3217*sb*sphi*st*vt)/100 + (3217*cb*cphi*sa*st*vt)/100)/vt;
dx6dot_dx5 = 0;
dx6dot_dx6 = (sb*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) - sb*vt*((15000*qbar*((15.0*P*(Cyp - delta_Cyp_lef*(lef/25.0 - 1.0)))/(vt*vt) + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(vt*vt)))/31847 + R*ca*cb - P*cb*sa) + ca*cb*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt) + cb*sa*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) - cb*sa*vt*(P*sb - Q*ca*cb + (84900*Q*qbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(31847.0*vt*vt)) - ca*cb*vt*(Q*cb*sa - R*sb + (84900*Q*qbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(31847*vt*vt)))/vt - (sb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + cb*sa*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt))/(vt*vt);
dx6dot_dx7 = -(Q*ca*ca*cb*cb*vt*vt - sb*vt*(P*ca*cb*vt + R*cb*sa*vt) + Q*cb*cb*sa*sa*vt*vt - ca*cb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + cb*sa*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt))/vt;
dx6dot_dx8 = -(sb*vt*(P*sa*sb*vt - R*ca*sb*vt) - cb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + sa*sb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + cb*sa*vt*(P*cb*vt + Q*ca*sb*vt) - ca*cb*vt*(R*cb*vt + Q*sa*sb*vt) + ca*sb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt))/vt;
dx6dot_dx9 = (sb*vt*((225000*qbar*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(31847*vt) + cb*sa*vt) - cb*sa*sb*vt*vt)/vt;
dx6dot_dx10 = (ca*cb*vt*((84900*qbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(31847*vt) - cb*sa*vt) + cb*sa*vt*((84900*qbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(31847*vt) + ca*cb*vt))/vt;
dx6dot_dx11 = (sb*vt*((225000*qbar*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(31847*vt) - ca*cb*vt) + ca*cb*sb*vt*vt)/vt;

// alpha_dot s relationship with other primary states
dx7dot_dx0 = 0;
dx7dot_dx1 = 0;
dx7dot_dx2 = 0;
dx7dot_dx3 = -(3217.0*ca*cb*ct*sphi*vt)/(100.0*(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt));
dx7dot_dx4 = ((3217.0*cb*ct*sa*vt)/100.0 - (3217.0*ca*cb*cphi*st*vt)/100)/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);
dx7dot_dx5 = 0;
dx7dot_dx6 = - (cb*sa*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt) - ca*cb*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*(P*sb - Q*ca*cb + (84900.0*Q*qbar*(Czq - delta_Cz_lef*(lef/25.0 - 1.0)))/(31847*vt*vt)) - cb*sa*vt*(Q*cb*sa - R*sb + (84900*Q*qbar*(Cxq - delta_Cxq_lef*(lef/25.0 - 1.0)))/(31847.0*vt*vt)))/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt) - ((2*vt*ca*ca*cb*cb + 2*vt*cb*cb*sa*sa)*(ca*cb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) - cb*sa*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)))/((ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt)*(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt));
dx7dot_dx7 = -(cb*sa*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt))/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);
dx7dot_dx8 = ((2*cb*sb*ca*ca*vt*vt + 2*cb*sb*sa*sa*vt*vt)*(ca*cb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) - cb*sa*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)))/((ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt)*(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt)) - (ca*sb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*(P*cb*vt + Q*ca*sb*vt) + cb*sa*vt*(R*cb*vt + Q*sa*sb*vt) - sa*sb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt))/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);
dx7dot_dx9 = -(ca*cb*sb*vt*vt)/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);
dx7dot_dx10 = (ca*cb*vt*((84900*qbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(31847*vt) + ca*cb*vt) - cb*sa*vt*((84900*qbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(31847*vt) - cb*sa*vt))/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);
dx7dot_dx11 = -(cb*sa*sb*vt*vt)/(ca*ca*cb*cb*vt*vt + cb*cb*sa*sa*vt*vt);

// beta_dot s relationship with other primary states
dx8dot_dx0 = 0;
dx8dot_dx1 = 0;
dx8dot_dx2 = 0;
dx8dot_dx3 = -(sb*((3217*cphi*ct*sb*vt)/100 - (3217*cb*ct*sa*sphi*vt)/100) - (3217*cphi*ct*vt)/100)/(cb*vt*vt);
dx8dot_dx4 = (sb*((3217*ca*cb*ct*vt)/100 + (3217*sb*sphi*st*vt)/100 + (3217*cb*cphi*sa*st*vt)/100) - (3217*sphi*st*vt)/100)/(cb*vt*vt);
dx8dot_dx5 = 0;
dx8dot_dx6 = (2*(sb*(sb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + cb*sa*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)) - vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt)))/(cb*vt*vt*vt) - (sb*(sb*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) - sb*vt*((15000*qbar*((15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(vt*vt)))/31847 + R*ca*cb - P*cb*sa) + ca*cb*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt) + cb*sa*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) - cb*sa*vt*(P*sb - Q*ca*cb + (84900*Q*qbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(31847*vt*vt)) - ca*cb*vt*(Q*cb*sa - R*sb + (84900*Q*qbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(31847*vt*vt))) - (3217*ct*sphi)/100 - (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 + vt*((15000*qbar*((15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(vt*vt)))/31847 + R*ca*cb - P*cb*sa) + R*ca*cb*vt - P*cb*sa*vt)/(cb*vt*vt);
dx8dot_dx7 = (vt*(P*ca*cb*vt + R*cb*sa*vt) + sb*(Q*ca*ca*cb*cb*vt*vt - sb*vt*(P*ca*cb*vt + R*cb*sa*vt) + Q*cb*cb*sa*sa*vt*vt - ca*cb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + cb*sa*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)))/(cb*vt*vt);
dx8dot_dx8 = - (vt*(P*sa*sb*vt - R*ca*sb*vt) - sb*(sb*vt*(P*sa*sb*vt - R*ca*sb*vt) - cb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + sa*sb*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + cb*sa*vt*(P*cb*vt + Q*ca*sb*vt) - ca*cb*vt*(R*cb*vt + Q*sa*sb*vt) + ca*sb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)) + cb*(sb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + cb*sa*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)))/(cb*vt*vt) - (sb*(sb*(sb*vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt) + cb*sa*vt*((3217*cphi*ct)/100 + (15000*qbar*(Cz - delta_Cz_lef*(lef/25 - 1) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1)))/(50*vt)))/31847 - P*sb*vt + Q*ca*cb*vt) + ca*cb*vt*((50*T)/31847 - (3217*st)/100 + (15000*qbar*(Cx - delta_Cx_lef*(lef/25 - 1) + (283*Q*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(50*vt)))/31847 + R*sb*vt - Q*cb*sa*vt)) - vt*((3217*ct*sphi)/100 + (15000*qbar*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/vt + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/vt))/31847 - R*ca*cb*vt + P*cb*sa*vt)))/(cb*cb*vt*vt);
dx8dot_dx9 = (vt*((225000*qbar*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(31847*vt) + cb*sa*vt) - sb*(sb*vt*((225000*qbar*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(31847*vt) + cb*sa*vt) - cb*sa*sb*vt*vt))/(cb*vt*vt);
dx8dot_dx10 = -(sb*(ca*cb*vt*((84900*qbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(31847*vt) - cb*sa*vt) + cb*sa*vt*((84900*qbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(31847*vt) + ca*cb*vt)))/(cb*vt*vt);
dx8dot_dx11 = (vt*((225000*qbar*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(31847*vt) - ca*cb*vt) - sb*(sb*vt*((225000*qbar*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(31847*vt) - ca*cb*vt) + ca*cb*sb*vt*vt))/(cb*vt*vt);

// P_dot s relationship with other primary states
dx9dot_dx0 = 0;
dx9dot_dx1 = 0;
dx9dot_dx2 = 0;
dx9dot_dx3 = 0;
dx9dot_dx4 = 0;
dx9dot_dx5 = 0;
dx9dot_dx6 = - (245500*qbar*((15*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(vt*vt) - (283*((15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(vt*vt))*(xcgr - 1/4))/750 + (15*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(vt*vt)))/16617591 - (15775000*qbar*((15*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(vt*vt)))/16617591;
dx9dot_dx7 = 0;
dx9dot_dx8 = (15775000*delta_Clbeta*qbar)/16617591 + (245500*delta_Cnbeta*qbar)/16617591;
dx9dot_dx9 = (1373327*Q)/49852773 + (245500*qbar*((15*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/vt - (283*(Cyp - delta_Cyp_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt)))/16617591 + (78875000*qbar*(Clp - delta_Clp_lef*(lef/25 - 1)))/(5539197*vt);
dx9dot_dx10 = (1373327*P)/49852773 - (38392577*R)/49852773;
dx9dot_dx11 = (245500*qbar*((15*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/vt - (283*(Cyr - delta_Cyr_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt)))/16617591 - (38392577*Q)/49852773 + (78875000*qbar*(Clr - delta_Clr_lef*(lef/25 - 1)))/(5539197*vt);

// Q_dot s relationship with other primary states
dx10dot_dx0 = 0;
dx10dot_dx1 = 0;
dx10dot_dx2 = 0;
dx10dot_dx3 = 0;
dx10dot_dx4 = 0;
dx10dot_dx5 = 0;
dx10dot_dx6 = -(1698*qbar*((283*Q*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(50*vt*vt) + (283*Q*(Czq - delta_Cz_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt*vt)))/27907;
dx10dot_dx7 = 0;
dx10dot_dx8 = 0;
dx10dot_dx9 = (26802*R)/27907 - (982*P)/27907;
dx10dot_dx10 = (1698*qbar*((283*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(50*vt) + (283*(Czq - delta_Cz_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt)))/27907;
dx10dot_dx11 = (26802*P)/27907 + (982*R)/27907;

// R_dot s relationship with other primary states
dx11dot_dx0 = 0;
dx11dot_dx1 = 0;
dx11dot_dx2 = 0;
dx11dot_dx3 = 0;
dx11dot_dx4 = 0;
dx11dot_dx5 = 0;
dx11dot_dx6 = - (245500*qbar*((15*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(vt*vt) - (283*((15*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(vt*vt))*(xcgr - 1/4))/750 + (15*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(vt*vt)))/16617591 - (15775000*qbar*((15*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(vt*vt) + (15*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(vt*vt)))/16617591;
dx11dot_dx7 = 0;
dx11dot_dx8 = (15775000*delta_Clbeta*qbar)/16617591 + (245500*delta_Cnbeta*qbar)/16617591;
dx11dot_dx9 = (1373327*Q)/49852773 + (245500*qbar*((15*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/vt - (283*(Cyp - delta_Cyp_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt)))/16617591 + (78875000*qbar*(Clp - delta_Clp_lef*(lef/25 - 1)))/(5539197*vt);
dx11dot_dx10 = (1373327*P)/49852773 - (38392577*R)/49852773;
dx11dot_dx11 = (245500*qbar*((15*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/vt - (283*(Cyr - delta_Cyr_lef*(lef/25 - 1))*(xcgr - 1/4))/(50*vt)))/16617591 - (38392577*Q)/49852773 + (78875000*qbar*(Clr - delta_Clr_lef*(lef/25 - 1)))/(5539197*vt);

// assign all this to jac
jac[0] = dx0dot_dx0;
jac[1] = dx0dot_dx1;
jac[2] = dx0dot_dx2;
jac[3] = dx0dot_dx3;
jac[4] = dx0dot_dx4;
jac[5] = dx0dot_dx5;
jac[6] = dx0dot_dx6;
jac[7] = dx0dot_dx7;
jac[8] = dx0dot_dx8;
jac[9] = dx0dot_dx9;
jac[10] = dx0dot_dx10;
jac[11] = dx0dot_dx11;

jac[12] = dx1dot_dx0;
jac[13] = dx1dot_dx1;
jac[14] = dx1dot_dx2;
jac[15] = dx1dot_dx3;
jac[16] = dx1dot_dx4;
jac[17] = dx1dot_dx5;
jac[18] = dx1dot_dx6;
jac[19] = dx1dot_dx7;
jac[20] = dx1dot_dx8;
jac[21] = dx1dot_dx9;
jac[22] = dx1dot_dx10;
jac[23] = dx1dot_dx11;

jac[24] = dx2dot_dx0;
jac[25] = dx2dot_dx1;
jac[26] = dx2dot_dx2;
jac[27] = dx2dot_dx3;
jac[28] = dx2dot_dx4;
jac[29] = dx2dot_dx5;
jac[30] = dx2dot_dx6;
jac[31] = dx2dot_dx7;
jac[32] = dx2dot_dx8;
jac[33] = dx2dot_dx9;
jac[34] = dx2dot_dx10;
jac[35] = dx2dot_dx11;

jac[36] = dx3dot_dx0;
jac[37] = dx3dot_dx1;
jac[38] = dx3dot_dx2;
jac[39] = dx3dot_dx3;
jac[40] = dx3dot_dx4;
jac[41] = dx3dot_dx5;
jac[42] = dx3dot_dx6;
jac[43] = dx3dot_dx7;
jac[44] = dx3dot_dx8;
jac[45] = dx3dot_dx9;
jac[46] = dx3dot_dx10;
jac[47] = dx3dot_dx11;

jac[48] = dx4dot_dx0;
jac[49] = dx4dot_dx1;
jac[50] = dx4dot_dx2;
jac[51] = dx4dot_dx3;
jac[52] = dx4dot_dx4;
jac[53] = dx4dot_dx5;
jac[54] = dx4dot_dx6;
jac[55] = dx4dot_dx7;
jac[56] = dx4dot_dx8;
jac[57] = dx4dot_dx9;
jac[58] = dx4dot_dx10;
jac[59] = dx4dot_dx11;

jac[60] = dx5dot_dx0;
jac[61] = dx5dot_dx1;
jac[62] = dx5dot_dx2;
jac[63] = dx5dot_dx3;
jac[64] = dx5dot_dx4;
jac[65] = dx5dot_dx5;
jac[66] = dx5dot_dx6;
jac[67] = dx5dot_dx7;
jac[68] = dx5dot_dx8;
jac[69] = dx5dot_dx9;
jac[70] = dx5dot_dx10;
jac[71] = dx5dot_dx11;

jac[72] = dx6dot_dx0;
jac[73] = dx6dot_dx1;
jac[74] = dx6dot_dx2;
jac[75] = dx6dot_dx3;
jac[76] = dx6dot_dx4;
jac[77] = dx6dot_dx5;
jac[78] = dx6dot_dx6;
jac[79] = dx6dot_dx7;
jac[80] = dx6dot_dx8;
jac[81] = dx6dot_dx9;
jac[82] = dx6dot_dx10;
jac[83] = dx6dot_dx11;

jac[84] = dx7dot_dx0;
jac[85] = dx7dot_dx1;
jac[86] = dx7dot_dx2;
jac[87] = dx7dot_dx3;
jac[88] = dx7dot_dx4;
jac[89] = dx7dot_dx5;
jac[90] = dx7dot_dx6;
jac[91] = dx7dot_dx7;
jac[92] = dx7dot_dx8;
jac[93] = dx7dot_dx9;
jac[94] = dx7dot_dx10;
jac[95] = dx7dot_dx11;

jac[96] = dx8dot_dx0;
jac[97] = dx8dot_dx1;
jac[98] = dx8dot_dx2;
jac[99] = dx8dot_dx3;
jac[100] = dx8dot_dx4;
jac[101] = dx8dot_dx5;
jac[102] = dx8dot_dx6;
jac[103] = dx8dot_dx7;
jac[104] = dx8dot_dx8;
jac[105] = dx8dot_dx9;
jac[106] = dx8dot_dx10;
jac[107] = dx8dot_dx11;

jac[108] = dx9dot_dx0;
jac[109] = dx9dot_dx1;
jac[110] = dx9dot_dx2;
jac[111] = dx9dot_dx3;
jac[112] = dx9dot_dx4;
jac[113] = dx9dot_dx5;
jac[114] = dx9dot_dx6;
jac[115] = dx9dot_dx7;
jac[116] = dx9dot_dx8;
jac[117] = dx9dot_dx9;
jac[118] = dx9dot_dx10;
jac[119] = dx9dot_dx11;

jac[120] = dx10dot_dx0;
jac[121] = dx10dot_dx1;
jac[122] = dx10dot_dx2;
jac[123] = dx10dot_dx3;
jac[124] = dx10dot_dx4;
jac[125] = dx10dot_dx5;
jac[126] = dx10dot_dx6;
jac[127] = dx10dot_dx7;
jac[128] = dx10dot_dx8;
jac[129] = dx10dot_dx9;
jac[130] = dx10dot_dx10;
jac[131] = dx10dot_dx11;

jac[132] = dx11dot_dx0;
jac[133] = dx11dot_dx1;
jac[134] = dx11dot_dx2;
jac[135] = dx11dot_dx3;
jac[136] = dx11dot_dx4;
jac[137] = dx11dot_dx5;
jac[138] = dx11dot_dx6;
jac[139] = dx11dot_dx7;
jac[140] = dx11dot_dx8;
jac[141] = dx11dot_dx9;
jac[142] = dx11dot_dx10;
jac[143] = dx11dot_dx11;



free(temp);

}; /*##### END of nlplant() ####*/