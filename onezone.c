#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <hdf5.h>
#include <grackle.h>

#define FAIL    0
#define SUCCESS 1
#define FLOAT double
#define Eint32 int

#define G00 1.7

#define MPROTON   1.67262171e-24
#define KBOLTZ    1.3806504e-16
#define GRAVITY   6.67428e-8
#define sec2Myr     3.1557e13

#define CM_PER_PC 3.085678e18
#define LSOLAR    3.85e33 // erg / s
#define hplanck      6.6260693e-27
#define ev2erg       1.60217653e-12
#define E_HI       (ev2erg * 13.6)
#define sigma_HI   6.30e-18
#define E_HeI      (ev2erg * 24.6)
#define sigma_HeI  7.42e-18
#define E_HeII     (ev2erg * 79.0)
#define sigma_HeII 1.575e-18
#define E_H2I      (ev2erg * 11.2)
#define sigma_H2I  3.83e-18

#define XH        0.76
#define Zsun      0.01295

#define tiny 1.0e-20

double Metallicity;
double Metallicity_loc;
double Metallicity_C13;
double Metallicity_C20;
double Metallicity_C25;
double Metallicity_C30;
double Metallicity_F13;
double Metallicity_F15;
double Metallicity_F50;
double Metallicity_F80;
double Metallicity_P170;
double Metallicity_P200;
double Metallicity_Y19;

//#define ISOTHERMAL
//#define GRAIN_GROWTH_TEST
#define ncell     1

enum ind_chem_species {
  i_elec, i_HI, i_HII, i_HeI, i_HeII, i_HeIII
, i_HM, i_H2I, i_H2II, i_DI, i_DII, i_HDI
, i_DM, i_HDII, i_HeHII
, i_CI, i_CII, i_CO, i_CO2, i_OI, i_OH, i_H2O, i_O2, i_SiI, i_SiOI, i_SiO2I
, i_CH, i_CH2, i_COII, i_OII, i_OHII, i_H2OII, i_H3OII, i_O2II
, i_Mg, i_Al, i_S, i_Fe
, i_SiM, i_FeM, i_Mg2SiO4, i_MgSiO3, i_Fe3O4, i_AC, i_SiO2D, i_MgO, i_FeS, i_Al2O3
, i_reforg, i_volorg, i_H2Oice
};

FILE *fdout;
char outfile[128];
void init_abundances(double *rho_gas, double *energy, double *gamma, double *rhoi);
void update     (double dt, double *density, double *pressure, double *energy
               , double *edot_adia, double edot_cool, double *gamma, double *rhoi
               , double *specific_heat_ratio, double *dt_col, int itchem);
void update_phys(double dt, double *density, double *pressure, double *energy
               , double *edot_adia, double edot_cool, double *gamma, double *rhoi
               , double *specific_heat_ratio, double *dt_col, int itchem);
double compute_collapse_time(double SHR, double rho_gas);
double timestep(int *itchem, double t_col, double t_chem, double t_cool, double dt);
double get_mu(double rhoi[]);
double get_ai(double temp, double rhoi[]);
void init_grackle(void);
int compute_abundances(int itime, double dtFixed, double rho, double *e, double *rhoi);
void output(int itime, double time, double rho, double energy, double rhoi[], double edot_adia);

code_units grackle_units;
grackle_field_data my_fields;
extern chemistry_data_storage grackle_rates;

double n0, T0;
int test, nZ, if0, iD0;

int main(int argc,char *argv[])
{

  int id = 0;
  int isp;
  int istep = 0;
  int itchem = 0;
  double rho_cgs, energy, energy_cool;
  double gamma, mu;
  double SHR, SHR_pre;
  double pressure, edot_adia, edot_cool;
  double dt, dt_current, t_chem, t_chem_i, t_col, t_cool, t_cool2;

  test = atoi(argv[1]);
  nZ   = atoi(argv[2]);

  if0  = 10;

  printf("[A] Grackle chemistry data structure.\n");
  chemistry_data *my_chemistry;
  my_chemistry = malloc(sizeof(chemistry_data));
  if (set_default_chemistry_parameters(my_chemistry) == FAIL) {
    printf("Error in grackle: set_default_chemistry_parameters\n");
    fflush(stdout);
    exit(1);
  }
  printf("[A] Done\n");

  grackle_data->primordial_chemistry = 4;

  n0 = 0.1;
  T0 = 300.0;

  Metallicity = pow(10.0, -nZ) * Zsun;
  Metallicity_loc = 1.0e-20 * Metallicity;
  Metallicity_C13 = 1.0e-20 * Metallicity;
  Metallicity_C20 = 1.0e-20 * Metallicity;
  Metallicity_C25 = 1.0e-20 * Metallicity;
  Metallicity_C30 = 0.5     * Metallicity;
  Metallicity_F13 = 1.0e-20 * Metallicity;
  Metallicity_F15 = 0.5     * Metallicity;
  Metallicity_F50 = 1.0e-20 * Metallicity;
  Metallicity_F80 = 1.0e-20 * Metallicity;
  Metallicity_P170= 1.0e-20 * Metallicity;
  Metallicity_P200= 1.0e-20 * Metallicity;
  Metallicity_Y19 = 1.0e-20 * Metallicity;
  grackle_data->SolarMetalFractionByMass = Zsun;

  switch(test) {
    case 0: /* Omukai (2005) */
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         = 0; // imabund loc
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 3; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 0; // igrgr
      grackle_data->dust_sublimation         = 1; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      break;
    case 1: /* Omukai et al. (2008) --- grackle still imcompatible with this test */
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         = 0; // imabund loc
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 3; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 0; // igrgr
      grackle_data->dust_sublimation         = 1; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      grackle_data->use_radiative_transfer = 1;
      /* TestA */
    //grackle_data->radiative_transfer_use_H2_shielding = 1;
    //grackle_data->H2_self_shielding                   = 0;
      /* TestB */
    //grackle_data->radiative_transfer_use_H2_shielding = 0;
    //grackle_data->H2_self_shielding                   = 1;
      /* TestC */
      grackle_data->radiative_transfer_use_H2_shielding = 0;
      grackle_data->H2_self_shielding                   = 3;
      break;
    case 2: /* Omukai et al. (2012) */
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         = 0; // imabund loc
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 3; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 0; // igrgr
      grackle_data->dust_sublimation         = 1; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      grackle_data->use_radiative_transfer = 1;
      /* TestA */
    //grackle_data->radiative_transfer_use_H2_shielding = 1;
    //grackle_data->H2_self_shielding                   = 0;
      /* TestB */
    //grackle_data->radiative_transfer_use_H2_shielding = 0;
    //grackle_data->H2_self_shielding                   = 1;
      /* TestC */
      grackle_data->radiative_transfer_use_H2_shielding = 0;
      grackle_data->H2_self_shielding                   = 3;
      T0 = 1000.0;
      iD0 = atoi(argv[3]);
      break;
    case 3: /* Marassi et al. (2014); Chiaki et al. (2015) */
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         = atoi(argv[3]); // imabund
      // 1: C13, 2: C20, 3: C25, 4: C30
      // 5: F13, 6: F15, 7: F50, 8: F80
      // 9: P170, 10: P200
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 2; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = atoi(argv[4]); // igrgr
      grackle_data->dust_sublimation         = 0; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      break;
    case 4: /* Park, Wise & Chiaki (2022) */
      grackle_data->metal_chemistry          = 0; // imchem
      grackle_data->metal_cooling            = 0; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         =11; // imabund Y19
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 1; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 0; // igrgr
      grackle_data->dust_sublimation         = 1; // idsub
      grackle_data->dust_temperature_multi   = 1; // itdmulti
      grackle_data->use_isrf_field               = 1; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      break;
    case 5: /* MultiMetals */
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 1; // immulti
      grackle_data->metal_abundances         = 0; // imabund
      grackle_data->use_dust_density_field   = 1; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 2; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 1; // igrgr
      grackle_data->dust_sublimation         = 0; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
      break;
    default:
      grackle_data->metal_chemistry          = 1; // imchem
      grackle_data->metal_cooling            = 1; // imcool
      grackle_data->multi_metals             = 0; // immulti
      grackle_data->metal_abundances         = 0; // imabund loc
      grackle_data->use_dust_density_field   = 0; // idustfield
      grackle_data->h2_on_dust               = 1; // idust
      grackle_data->dust_species             = 0; // idspecies
      grackle_data->dust_chemistry           = 0; // idustall
      grackle_data->grain_growth             = 0; // igrgr
      grackle_data->dust_sublimation         = 0; // idsub
      grackle_data->dust_temperature_multi   = 0; // itdmulti
      grackle_data->use_isrf_field               = 0; // iisrffield
      grackle_data->interstellar_radiation_field = 0; // isrf
  }

  grackle_data->HydrogenFractionByMass = XH;
  grackle_data->grackle_data_file = "CloudyData_noUVB.h5";


  /* Others */
  grackle_data->use_grackle                    = 1;
#ifdef ISOTHERMAL
  grackle_data->with_radiative_cooling         = 0;
#else
  grackle_data->with_radiative_cooling         = 1;
#endif
//sprintf(grackle_data->grackle_data_file, "./GRACKLE_DATA_FILE");
  grackle_data->UVbackground                   = 0;
  grackle_data->Compton_xray_heating           = 0;
  grackle_data->LWbackground_intensity         = 0.0;
  grackle_data->LWbackground_sawtooth_suppression = 0;
  grackle_data->Gamma                          = 5.0/3.0;
//grackle_data->primordial_chemistry           = 19;
//grackle_data->metal_cooling                  = 0;
//grackle_data->h2_on_dust                     = 1;
  grackle_data->cmb_temperature_floor          = 1;
  grackle_data->three_body_rate                = 1; // PSS83
//grackle_data->three_body_rate                = 4; // GA08
//grackle_data->three_body_rate                = 5; // Forrey
  grackle_data->cie_cooling                    = 1;
  grackle_data->h2_optical_depth_approximation = 1;
  grackle_data->photoelectric_heating          = 0;
  grackle_data->photoelectric_heating_rate     = 0.0;
//grackle_data->NumberOfTemperatureBins        = 901;
  grackle_data->CaseBRecombination             = 1;
//grackle_data->TemperatureStart               = 1.0;
//grackle_data->TemperatureEnd                 = 1.0e9;
//grackle_data->NumberOfDustTemperatureBins    = 901;
//grackle_data->DustTemperatureStart           = 1.0;
//grackle_data->DustTemperatureEnd             = 1.0e4;
//grackle_data->HydrogenFractionByMass         = XH;
//grackle_data->DeuteriumToHydrogenRatio       = 3.0e-5*2.0;
//grackle_data->SolarMetalFractionByMass       = Zsun;
//grackle_data->MetalFractionByMass            = 1.0e-4;
//grackle_data->metal_chemistry                = 1;
//grackle_data->grain_growth                   = 0;
///************************ RECOMMEND ***********************/
    grackle_data->use_palla_salpeter_stahler_1983        = 1;
    grackle_data->use_stancil_lepp_dalgarno_1998         = 1;
    grackle_data->use_omukai_gas_grain                   = 1;
    grackle_data->use_uniform_grain_dist_gamma_isrf      = 1;
///**********************************************************/

  // Initialize units structure.
  grackle_units.a_units              = 1.0;
  if(test == 0)
      grackle_units.a_value              = 1.0 / grackle_units.a_units;
  else
      grackle_units.a_value              = 0.1 / grackle_units.a_units;
  grackle_units.comoving_coordinates = 0;
  grackle_units.density_units        = 1.0e0 * MPROTON / grackle_data->HydrogenFractionByMass;
  grackle_units.length_units         = sqrt(M_PI * 5.0/3.0 * KBOLTZ * 300.0 / 1.23 / MPROTON / GRAVITY / grackle_units.density_units);
  grackle_units.time_units           = sqrt(3.0 * M_PI / 32.0 / GRAVITY / grackle_units.density_units);
  grackle_units.velocity_units       = grackle_units.length_units / grackle_units.time_units;

  grackle_verbose = 1;

  printf("[B] Initialize chemistry structure.\n");
  if (initialize_chemistry_data(&grackle_units) == FAIL) {
    printf("Error in Grackle initialize_chemistry_data.\n");
    exit(1);
  }
  printf("[B] Done\n");

  init_grackle();

  int n_sp = 51;
  double rhoi[n_sp], rhoi_pre[n_sp];

  printf("[C] init_abundances\n");
  init_abundances(&rho_cgs, &energy, &gamma, rhoi);
  printf("[C] done\n");

  switch(test) {
    case 0: /* Omukai (2005) */
      sprintf(outfile, "O05/output_Z-%d", nZ);
      break;
    case 1: /* Omukai et al. (2008) --- grackle still imcompatible with this test */
      sprintf(outfile, "O08/output_Z-%d", nZ);
      break;
    case 2: /* Omukai et al. (2012) */
      sprintf(outfile, "O12/output_Z-%d_D-%d", nZ, iD0);
      break;
    case 3: /* Marassi et al. (2014); Chiaki et al. (2015) */
      if(atoi(argv[3]) >= 5 && atoi(argv[3]) <= 8)
        sprintf(outfile, "M14/output_SN%02d_Z-%d_gg%d", grackle_data->metal_abundances, nZ, grackle_data->grain_growth);
      else
        sprintf(outfile, "C15/output_SN%02d_Z-%d_gg%d", grackle_data->metal_abundances, nZ, grackle_data->grain_growth);
      break;
    case 4: /* Park, Wise & Chiaki (2022) */
      sprintf(outfile, "P22/output_Z-%d", nZ);
      break;
    case 5: /* MultiMetals */
      sprintf(outfile, "multimetal/output_Z-%d", nZ);
      break;
    default:
      sprintf(outfile, "default/output_Z-%d", nZ);
  }

  fdout=fopen(outfile, "w");

  dt  = 0.1;
  SHR = 5.0/3.0;
  pressure = (gamma-1.0)*energy;
  t_col    = sqrt(3.0*M_PI/32.0/GRAVITY/rho_cgs);
  edot_adia = pressure/t_col;
  itchem = 0;

//printf("[%s] %3d i=%9d z=%13.5e dt=%13.5e d=%13.5e u=%13.5e g=%13.5e dv=%13.5e\n"
//  , __FUNCTION__, istep, id, redshift, dt, rho_cgs, energy, gamma, dvdr);
//for(isp=0; isp<n_sp; isp++) printf("%9.2e ", rhoi[isp]); printf("\n");

  double cool_rate;

  int itime;
  double time;
  double t0, t1, tc0, tc1, t_comp;

  t0 = ((double) clock()) / CLOCKS_PER_SEC;
  t_comp = 0.0;

  time = 0.0;
  itime = 0;
  do {

    for(isp=0; isp<n_sp; isp++) rhoi_pre[isp] = rhoi[isp];
    energy_cool = energy;
 
    tc0 = ((double) clock());
    if(compute_abundances(itime, dt, rho_cgs, &energy_cool, rhoi) == FAIL) {
      exit(1);
    }
    tc1 = ((double) clock());
    t_comp += tc1 - tc0;
    if(ncell > 1) {
      printf("%13.5e %23.15e\n", XH*rho_cgs/1.6726e-24, t_comp / CLOCKS_PER_SEC);
      fflush(stdout);
    }
//printf("A %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

//  double t_col;
//  t_col = 4.7*sqrt(1.0/24.0/M_PI/GRAVITY/rho_cgs);
//  edot_cool = - 0.3* (5.0/3.0-1.0) * energy/t_col;
//  energy_cool = energy + edot_cool * dt;

    t_chem = 1.0e0;
    for(isp=0; isp<n_sp; isp++) {
      t_chem_i = 0.5*(rhoi[isp] + rhoi_pre[isp]) / fabs(rhoi[isp] - rhoi_pre[isp]) * dt;
      if(t_chem_i && t_chem_i < t_chem) t_chem = t_chem_i;
    }
//printf("B %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

#if defined(ISOTHERMAL) || defined(GRAIN_GROWTH_TEST)
    edot_cool = 1.0e-100;
#else
    edot_cool = (energy_cool - energy) / dt;
#endif
//printf("C %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

//  t_cool = 1.0e100; ///// TEMPO /////
    t_cool = 0.5*(energy_cool + energy) / edot_cool;
//  t_cool2 = compute_cooling_time(1.0/(1.0+redshift), dt, rho_cgs, &energy, rhoi);
//printf("D %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

//  update     (dt, &rho_cgs, &pressure, &energy, &edot_adia, edot_cool, &gamma, rhoi, &SHR, &t_col, itchem);
    update_phys(dt, &rho_cgs, &pressure, &energy, &edot_adia, edot_cool, &gamma, rhoi, &SHR, &t_col, itchem);
    /* In update_phys(), energy is updated from the one before radiative cooling. */
//printf("E %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

//  printf("UPDATE %13.5e %13.5e %13.5e %13.5e %13.5e\n", rho_cgs*XH/1.67e-24, edot_adia, edot_cool,dt, energy);
//  printf("UPDATE %13.5e %13.5e\n", get_mu(rhoi), gamma);

    dt_current = dt;
    dt = timestep(&itchem, t_col, t_chem, t_cool, dt_current);
//printf("F %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

    output(itime, time, rho_cgs, energy, rhoi, edot_adia);
//printf("G %13.5e %13.5e \n", rhoi[i_HI]/rho_cgs, rhoi[i_OH]/rho_cgs);

    if(XH * rho_cgs / MPROTON > 1.0e3) {
      time += dt;
//    printf("%13.5e %13.5e\n", time / sec2Myr, XH * rho_cgs / MPROTON);
    }
    itime++;

  } while(
#ifdef GRAIN_GROWTH_TEST
      time < 1.0e11
#else
      rho_cgs < 1.0e16*MPROTON/grackle_data->HydrogenFractionByMass
#endif
  );

//t1 = ((double) clock()) / CLOCKS_PER_SEC;
//t_comp = t1 - t0;

  fclose(fdout);

  if(ncell == 1) {
    t_comp /= CLOCKS_PER_SEC;
    printf("COMP TIME %13.7f sec %9d loops %13.7f sec/loop\n", t_comp, itime, t_comp/(double)itime);
    printf("Elapsed time %13.7f Myr\n", time / sec2Myr);
  }

  return 0;
}


void init_abundances(double *rho_gas, double *energy, double *gamma, double *rhoi)
{
  int isp;
  int n_sp = 6 + grackle_data->primordial_chemistry * 3;
  n_sp = 51;
  double mu;
  double nHmH;
  double H_frac, D_frac, He_frac;
  double C_frac, O_frac, Mg_frac, Al_frac, Si_frac, S_frac, Fe_frac;
  double metal_density, dust_density;
  int    iSN0;
  double metal_loc, metal_C13, metal_C20, metal_C25, metal_C30
       , metal_F13, metal_F15, metal_F50, metal_F80
       , metal_P170, metal_P200, metal_Y19;
  double SiM_frac     , FeM_frac    , Mg2SiO4_frac, MgSiO3_frac , Fe3O4_frac  
       , AC_frac      , SiO2D_frac  , MgO_frac    , FeS_frac    , Al2O3_frac
       , reforg_frac  , volorg_frac , H2Oice_frac ;

#ifdef GRAIN_GROWTH_TEST
  n0 = 1.0e12;
#else
  n0 = 0.1;
//n0 = 1.0e5;
//n0 = 1.0e12;
#endif
#ifdef ISOTHERMAL
  T0 = 2e3;
//#else
//#ifdef RADIATION
//  T0 = 1000.0;
//#else
//  T0 = 300.0;
//#endif
#endif
//n0 = 1.e6;
//T0 = 1.e4;

  for(isp=0; isp<n_sp; isp++) rhoi[isp] = 0.0;

   H_frac = (1.0 - grackle_data->DeuteriumToHydrogenRatio) * grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
// H_frac = grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
   D_frac = grackle_data->DeuteriumToHydrogenRatio * grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
  He_frac = (1.0 - grackle_data->HydrogenFractionByMass) * (1.0 - Metallicity);

  if(grackle_data->multi_metals == 1) {
    metal_loc     = Metallicity_loc;
    metal_C13     = Metallicity_C13;
    metal_C20     = Metallicity_C20;
    metal_C25     = Metallicity_C25;
    metal_C30     = Metallicity_C30;
    metal_F13     = Metallicity_F13;
    metal_F15     = Metallicity_F15;
    metal_F50     = Metallicity_F50;
    metal_F80     = Metallicity_F80;
    metal_P170    = Metallicity_P170;
    metal_P200    = Metallicity_P200;
    metal_Y19     = Metallicity_Y19;
    metal_density = Metallicity_loc
                  + Metallicity_C13
                  + Metallicity_C20
                  + Metallicity_C25
                  + Metallicity_C30
                  + Metallicity_F13
                  + Metallicity_F15
                  + Metallicity_F50
                  + Metallicity_F80
                  + Metallicity_P170
                  + Metallicity_P200
                  + Metallicity_Y19;
  }
  dust_density = 0.0;

  if(grackle_data->multi_metals == 0) {

     iSN0 = grackle_data->metal_abundances;

//  if(grackle_data->grain_growth) {
     C_frac  = grackle_data->SN0_fC [iSN0] * Metallicity;
     O_frac  = grackle_data->SN0_fO [iSN0] * Metallicity;
    Mg_frac  = grackle_data->SN0_fMg[iSN0] * Metallicity;
    Al_frac  = grackle_data->SN0_fAl[iSN0] * Metallicity;
    Si_frac  = grackle_data->SN0_fSi[iSN0] * Metallicity;
     S_frac  = grackle_data->SN0_fS [iSN0] * Metallicity;
    Fe_frac  = grackle_data->SN0_fFe[iSN0] * Metallicity;
  
        SiM_frac = grackle_data->SN0_fSiM    [iSN0] * Metallicity;
        FeM_frac = grackle_data->SN0_fFeM    [iSN0] * Metallicity;
    Mg2SiO4_frac = grackle_data->SN0_fMg2SiO4[iSN0] * Metallicity;
     MgSiO3_frac = grackle_data->SN0_fMgSiO3 [iSN0] * Metallicity;
      Fe3O4_frac = grackle_data->SN0_fFe3O4  [iSN0] * Metallicity;
         AC_frac = grackle_data->SN0_fAC     [iSN0] * Metallicity;
      SiO2D_frac = grackle_data->SN0_fSiO2D  [iSN0] * Metallicity;
        MgO_frac = grackle_data->SN0_fMgO    [iSN0] * Metallicity;
        FeS_frac = grackle_data->SN0_fFeS    [iSN0] * Metallicity;
      Al2O3_frac = grackle_data->SN0_fAl2O3  [iSN0] * Metallicity;
     reforg_frac = grackle_data->SN0_freforg [iSN0] * Metallicity;
     volorg_frac = grackle_data->SN0_fvolorg [iSN0] * Metallicity;
     H2Oice_frac = grackle_data->SN0_fH2Oice [iSN0] * Metallicity;

//  } else {
//   C_frac  = grackle_data->SN0_XC [iSN0] * Metallicity;
//   O_frac  = grackle_data->SN0_XO [iSN0] * Metallicity;
//  Mg_frac  = grackle_data->SN0_XMg[iSN0] * Metallicity;
//  Al_frac  = grackle_data->SN0_XAl[iSN0] * Metallicity;
//  Si_frac  = grackle_data->SN0_XSi[iSN0] * Metallicity;
//   S_frac  = grackle_data->SN0_XS [iSN0] * Metallicity;
//  Fe_frac  = grackle_data->SN0_XFe[iSN0] * Metallicity;
//  }

  } else {

     C_frac  = 0.0;
     O_frac  = 0.0;
    Mg_frac  = 0.0;
    Al_frac  = 0.0;
    Si_frac  = 0.0;
     S_frac  = 0.0;
    Fe_frac  = 0.0;
  
     C_frac += grackle_data->SN0_fC [0] * metal_loc;
     O_frac += grackle_data->SN0_fO [0] * metal_loc;
    Mg_frac += grackle_data->SN0_fMg[0] * metal_loc;
    Al_frac += grackle_data->SN0_fAl[0] * metal_loc;
    Si_frac += grackle_data->SN0_fSi[0] * metal_loc;
     S_frac += grackle_data->SN0_fS [0] * metal_loc;
    Fe_frac += grackle_data->SN0_fFe[0] * metal_loc;
  
     C_frac += grackle_data->SN0_fC [1] * metal_C13;
     O_frac += grackle_data->SN0_fO [1] * metal_C13;
    Mg_frac += grackle_data->SN0_fMg[1] * metal_C13;
    Al_frac += grackle_data->SN0_fAl[1] * metal_C13;
    Si_frac += grackle_data->SN0_fSi[1] * metal_C13;
     S_frac += grackle_data->SN0_fS [1] * metal_C13;
    Fe_frac += grackle_data->SN0_fFe[1] * metal_C13;
  
     C_frac += grackle_data->SN0_fC [2] * metal_C20;
     O_frac += grackle_data->SN0_fO [2] * metal_C20;
    Mg_frac += grackle_data->SN0_fMg[2] * metal_C20;
    Al_frac += grackle_data->SN0_fAl[2] * metal_C20;
    Si_frac += grackle_data->SN0_fSi[2] * metal_C20;
     S_frac += grackle_data->SN0_fS [2] * metal_C20;
    Fe_frac += grackle_data->SN0_fFe[2] * metal_C20;
  
     C_frac += grackle_data->SN0_fC [3] * metal_C25;
     O_frac += grackle_data->SN0_fO [3] * metal_C25;
    Mg_frac += grackle_data->SN0_fMg[3] * metal_C25;
    Al_frac += grackle_data->SN0_fAl[3] * metal_C25;
    Si_frac += grackle_data->SN0_fSi[3] * metal_C25;
     S_frac += grackle_data->SN0_fS [3] * metal_C25;
    Fe_frac += grackle_data->SN0_fFe[3] * metal_C25;
  
     C_frac += grackle_data->SN0_fC [4] * metal_C30;
     O_frac += grackle_data->SN0_fO [4] * metal_C30;
    Mg_frac += grackle_data->SN0_fMg[4] * metal_C30;
    Al_frac += grackle_data->SN0_fAl[4] * metal_C30;
    Si_frac += grackle_data->SN0_fSi[4] * metal_C30;
     S_frac += grackle_data->SN0_fS [4] * metal_C30;
    Fe_frac += grackle_data->SN0_fFe[4] * metal_C30;
  
     C_frac += grackle_data->SN0_fC [5] * metal_F13;
     O_frac += grackle_data->SN0_fO [5] * metal_F13;
    Mg_frac += grackle_data->SN0_fMg[5] * metal_F13;
    Al_frac += grackle_data->SN0_fAl[5] * metal_F13;
    Si_frac += grackle_data->SN0_fSi[5] * metal_F13;
     S_frac += grackle_data->SN0_fS [5] * metal_F13;
    Fe_frac += grackle_data->SN0_fFe[5] * metal_F13;
  
     C_frac += grackle_data->SN0_fC [6] * metal_F15;
     O_frac += grackle_data->SN0_fO [6] * metal_F15;
    Mg_frac += grackle_data->SN0_fMg[6] * metal_F15;
    Al_frac += grackle_data->SN0_fAl[6] * metal_F15;
    Si_frac += grackle_data->SN0_fSi[6] * metal_F15;
     S_frac += grackle_data->SN0_fS [6] * metal_F15;
    Fe_frac += grackle_data->SN0_fFe[6] * metal_F15;
  
     C_frac += grackle_data->SN0_fC [7] * metal_F50;
     O_frac += grackle_data->SN0_fO [7] * metal_F50;
    Mg_frac += grackle_data->SN0_fMg[7] * metal_F50;
    Al_frac += grackle_data->SN0_fAl[7] * metal_F50;
    Si_frac += grackle_data->SN0_fSi[7] * metal_F50;
     S_frac += grackle_data->SN0_fS [7] * metal_F50;
    Fe_frac += grackle_data->SN0_fFe[7] * metal_F50;
  
     C_frac += grackle_data->SN0_fC [8] * metal_F80;
     O_frac += grackle_data->SN0_fO [8] * metal_F80;
    Mg_frac += grackle_data->SN0_fMg[8] * metal_F80;
    Al_frac += grackle_data->SN0_fAl[8] * metal_F80;
    Si_frac += grackle_data->SN0_fSi[8] * metal_F80;
     S_frac += grackle_data->SN0_fS [8] * metal_F80;
    Fe_frac += grackle_data->SN0_fFe[8] * metal_F80;
  
     C_frac += grackle_data->SN0_fC [9] * metal_P170;
     O_frac += grackle_data->SN0_fO [9] * metal_P170;
    Mg_frac += grackle_data->SN0_fMg[9] * metal_P170;
    Al_frac += grackle_data->SN0_fAl[9] * metal_P170;
    Si_frac += grackle_data->SN0_fSi[9] * metal_P170;
     S_frac += grackle_data->SN0_fS [9] * metal_P170;
    Fe_frac += grackle_data->SN0_fFe[9] * metal_P170;
  
     C_frac += grackle_data->SN0_fC [10]* metal_P200;
     O_frac += grackle_data->SN0_fO [10]* metal_P200;
    Mg_frac += grackle_data->SN0_fMg[10]* metal_P200;
    Al_frac += grackle_data->SN0_fAl[10]* metal_P200;
    Si_frac += grackle_data->SN0_fSi[10]* metal_P200;
     S_frac += grackle_data->SN0_fS [10]* metal_P200;
    Fe_frac += grackle_data->SN0_fFe[10]* metal_P200;
  
     C_frac += grackle_data->SN0_fC [11]* metal_Y19;
     O_frac += grackle_data->SN0_fO [11]* metal_Y19;
    Mg_frac += grackle_data->SN0_fMg[11]* metal_Y19;
    Al_frac += grackle_data->SN0_fAl[11]* metal_Y19;
    Si_frac += grackle_data->SN0_fSi[11]* metal_Y19;
     S_frac += grackle_data->SN0_fS [11]* metal_Y19;
    Fe_frac += grackle_data->SN0_fFe[11]* metal_Y19;
  
        SiM_frac  = 0.0;
        FeM_frac  = 0.0;
    Mg2SiO4_frac  = 0.0;
     MgSiO3_frac  = 0.0;
      Fe3O4_frac  = 0.0;
         AC_frac  = 0.0;
      SiO2D_frac  = 0.0;
        MgO_frac  = 0.0;
        FeS_frac  = 0.0;
      Al2O3_frac  = 0.0;
     reforg_frac  = 0.0;
     volorg_frac  = 0.0;
     H2Oice_frac  = 0.0;
  
        FeM_frac += grackle_data->SN0_fFeM    [0] * Metallicity_loc;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[0] * Metallicity_loc;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [0] * Metallicity_loc;
        FeS_frac += grackle_data->SN0_fFeS    [0] * Metallicity_loc;
     reforg_frac += grackle_data->SN0_freforg [0] * Metallicity_loc;
     volorg_frac += grackle_data->SN0_fvolorg [0] * Metallicity_loc;
     H2Oice_frac += grackle_data->SN0_fH2Oice [0] * Metallicity_loc;
  
        SiM_frac += grackle_data->SN0_fSiM    [1] * Metallicity_C13;
        FeM_frac += grackle_data->SN0_fFeM    [1] * Metallicity_C13;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[1] * Metallicity_C13;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [1] * Metallicity_C13;
         AC_frac += grackle_data->SN0_fAC     [1] * Metallicity_C13;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [1] * Metallicity_C13;
        MgO_frac += grackle_data->SN0_fMgO    [1] * Metallicity_C13;
        FeS_frac += grackle_data->SN0_fFeS    [1] * Metallicity_C13;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [1] * Metallicity_C13;
  
        SiM_frac += grackle_data->SN0_fSiM    [2] * Metallicity_C20;
        FeM_frac += grackle_data->SN0_fFeM    [2] * Metallicity_C20;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[2] * Metallicity_C20;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [2] * Metallicity_C20;
         AC_frac += grackle_data->SN0_fAC     [2] * Metallicity_C20;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [2] * Metallicity_C20;
        MgO_frac += grackle_data->SN0_fMgO    [2] * Metallicity_C20;
        FeS_frac += grackle_data->SN0_fFeS    [2] * Metallicity_C20;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [2] * Metallicity_C20;
  
        SiM_frac += grackle_data->SN0_fSiM    [3] * Metallicity_C25;
        FeM_frac += grackle_data->SN0_fFeM    [3] * Metallicity_C25;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[3] * Metallicity_C25;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [3] * Metallicity_C25;
         AC_frac += grackle_data->SN0_fAC     [3] * Metallicity_C25;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [3] * Metallicity_C25;
        MgO_frac += grackle_data->SN0_fMgO    [3] * Metallicity_C25;
        FeS_frac += grackle_data->SN0_fFeS    [3] * Metallicity_C25;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [3] * Metallicity_C25;
  
        SiM_frac += grackle_data->SN0_fSiM    [4] * Metallicity_C30;
        FeM_frac += grackle_data->SN0_fFeM    [4] * Metallicity_C30;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[4] * Metallicity_C30;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [4] * Metallicity_C30;
         AC_frac += grackle_data->SN0_fAC     [4] * Metallicity_C30;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [4] * Metallicity_C30;
        MgO_frac += grackle_data->SN0_fMgO    [4] * Metallicity_C30;
        FeS_frac += grackle_data->SN0_fFeS    [4] * Metallicity_C30;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [4] * Metallicity_C30;
  
        FeM_frac += grackle_data->SN0_fFeM    [5] * Metallicity_F13;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[5] * Metallicity_F13;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [5] * Metallicity_F13;
      Fe3O4_frac += grackle_data->SN0_fFe3O4  [5] * Metallicity_F13;
         AC_frac += grackle_data->SN0_fAC     [5] * Metallicity_F13;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [5] * Metallicity_F13;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [5] * Metallicity_F13;
  
        FeM_frac += grackle_data->SN0_fFeM    [6] * Metallicity_F15;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[6] * Metallicity_F15;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [6] * Metallicity_F15;
      Fe3O4_frac += grackle_data->SN0_fFe3O4  [6] * Metallicity_F15;
         AC_frac += grackle_data->SN0_fAC     [6] * Metallicity_F15;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [6] * Metallicity_F15;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [6] * Metallicity_F15;
  
        FeM_frac += grackle_data->SN0_fFeM    [7] * Metallicity_F50;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[7] * Metallicity_F50;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [7] * Metallicity_F50;
      Fe3O4_frac += grackle_data->SN0_fFe3O4  [7] * Metallicity_F50;
         AC_frac += grackle_data->SN0_fAC     [7] * Metallicity_F50;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [7] * Metallicity_F50;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [7] * Metallicity_F50;
  
        FeM_frac += grackle_data->SN0_fFeM    [8] * Metallicity_F80;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[8] * Metallicity_F80;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [8] * Metallicity_F80;
      Fe3O4_frac += grackle_data->SN0_fFe3O4  [8] * Metallicity_F80;
         AC_frac += grackle_data->SN0_fAC     [8] * Metallicity_F80;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [8] * Metallicity_F80;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [8] * Metallicity_F80;

        SiM_frac += grackle_data->SN0_fSiM    [9] * Metallicity_P170;
        FeM_frac += grackle_data->SN0_fFeM    [9] * Metallicity_P170;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[9] * Metallicity_P170;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [9] * Metallicity_P170;
         AC_frac += grackle_data->SN0_fAC     [9] * Metallicity_P170;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [9] * Metallicity_P170;
        MgO_frac += grackle_data->SN0_fMgO    [9] * Metallicity_P170;
        FeS_frac += grackle_data->SN0_fFeS    [9] * Metallicity_P170;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [9] * Metallicity_P170;
  
        SiM_frac += grackle_data->SN0_fSiM    [10]* Metallicity_P200;
        FeM_frac += grackle_data->SN0_fFeM    [10]* Metallicity_P200;
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[10]* Metallicity_P200;
     MgSiO3_frac += grackle_data->SN0_fMgSiO3 [10]* Metallicity_P200;
         AC_frac += grackle_data->SN0_fAC     [10]* Metallicity_P200;
      SiO2D_frac += grackle_data->SN0_fSiO2D  [10]* Metallicity_P200;
        MgO_frac += grackle_data->SN0_fMgO    [10]* Metallicity_P200;
        FeS_frac += grackle_data->SN0_fFeS    [10]* Metallicity_P200;
      Al2O3_frac += grackle_data->SN0_fAl2O3  [10]* Metallicity_P200;
                                                 
    Mg2SiO4_frac += grackle_data->SN0_fMg2SiO4[11]* Metallicity_Y19;
         AC_frac += grackle_data->SN0_fAC     [11]* Metallicity_Y19;
  

  }

  printf("   C_frac %13.7f\n", 12 + log10( C_frac / XH / 12.0) - 8.43);
  printf("   O_frac %13.7f\n", 12 + log10( O_frac / XH / 16.0) - 8.69);
  printf("  Mg_frac %13.7f\n", 12 + log10(Mg_frac / XH / 24.0) - 7.60);
  printf("  Al_frac %13.7f\n", 12 + log10(Al_frac / XH / 27.0) - 6.45);
  printf("  Si_frac %13.7f\n", 12 + log10(Si_frac / XH / 28.0) - 7.51);
  printf("   S_frac %13.7f\n", 12 + log10( S_frac / XH / 32.0) - 7.12);
  printf("  Fe_frac %13.7f\n", 12 + log10(Fe_frac / XH / 56.0) - 7.50);

  rhoi[i_CI     ] = tiny; 
  rhoi[i_CII    ] = C_frac;
  rhoi[i_CO     ] = tiny; 
  rhoi[i_CO2    ] = tiny; 
  rhoi[i_OI     ] = O_frac;
  rhoi[i_OH     ] = tiny; 
  rhoi[i_H2O    ] = tiny; 
  rhoi[i_O2     ] = tiny; 
  rhoi[i_SiI    ] = Si_frac;
  rhoi[i_SiOI   ] = tiny; 
  rhoi[i_SiO2I  ] = tiny; 
  rhoi[i_CH     ] = tiny; 
  rhoi[i_CH2    ] = tiny; 
  rhoi[i_COII   ] = tiny; 
  rhoi[i_OII    ] = tiny; 
  rhoi[i_OHII   ] = tiny; 
  rhoi[i_H2OII  ] = tiny; 
  rhoi[i_H3OII  ] = tiny; 
  rhoi[i_O2II   ] = tiny; 
  rhoi[i_Mg     ] = Mg_frac;
  rhoi[i_Al     ] = Al_frac;
  rhoi[i_S      ] = (isnan(S_frac) ? tiny : S_frac);
  rhoi[i_Fe     ] = Fe_frac;

  rhoi[i_SiM     ] =      SiM_frac;
  rhoi[i_FeM     ] =      FeM_frac;
  rhoi[i_Mg2SiO4 ] =  Mg2SiO4_frac;
  rhoi[i_MgSiO3  ] =   MgSiO3_frac;
  rhoi[i_Fe3O4   ] =    Fe3O4_frac;
  rhoi[i_AC      ] =       AC_frac;
  rhoi[i_SiO2D   ] =    SiO2D_frac;
  rhoi[i_MgO     ] =      MgO_frac;
  rhoi[i_FeS     ] =      FeS_frac;
  rhoi[i_Al2O3   ] =    Al2O3_frac;
  rhoi[i_reforg  ] =   reforg_frac;
  rhoi[i_volorg  ] =   volorg_frac;
  rhoi[i_H2Oice  ] =   H2Oice_frac;

  rhoi[i_HeHII  ] = He_frac * 1.0e-20;
  rhoi[i_HeII   ] = He_frac * 1.0e-20;
  rhoi[i_HeIII  ] = He_frac * 1.0e-20;
  rhoi[i_HeI    ] = He_frac
                  - rhoi[i_HeII] 
                  - rhoi[i_HeIII];
  if(grackle_data->primordial_chemistry > 3)
    rhoi[i_HeI] = rhoi[i_HeI]
                  - rhoi[i_HeHII] * 4.0/5.0;

  rhoi[i_DII    ] = D_frac * 1.0e-4;
  rhoi[i_HDI    ] = D_frac * 1.0e-6;
  rhoi[i_DM     ] = D_frac * 1.0e-20;
  rhoi[i_HDII   ] = D_frac * 1.0e-20;
  rhoi[i_DI     ] = D_frac
                  - rhoi[i_DII]
                  - rhoi[i_HDI] * 2.0/3.0;
  if(grackle_data->primordial_chemistry > 3)
    rhoi[i_DI] = rhoi[i_DI]
                  - rhoi[i_DM]
                  - rhoi[i_HDII] * 2.0/3.0;

  rhoi[i_HII    ] = H_frac * 1.0e-4;
  rhoi[i_H2I    ] = H_frac * 1.0e-6;
  rhoi[i_HM     ] = H_frac * 1.0e-20;
  rhoi[i_H2II   ] = H_frac * 1.0e-20;
  rhoi[i_HI     ] = H_frac - rhoi[i_HII];
  if(grackle_data->primordial_chemistry > 1)
    rhoi[i_HI] = rhoi[i_HI]
                  - rhoi[i_H2I]
                  - rhoi[i_HM]
                  - rhoi[i_H2II];
  if(grackle_data->primordial_chemistry > 2)
    rhoi[i_HI] = rhoi[i_HI]
                  - rhoi[i_HDI]   * 1.0/3.0;
  if(grackle_data->primordial_chemistry > 3)
       rhoi[i_HI] = rhoi[i_HI]
                  - rhoi[i_HDII]  * 1.0/3.0
                  - rhoi[i_HeHII] * 1.0/5.0;
  if(grackle_data->metal_chemistry)
    if(Metallicity > 1.0e-9) 
      rhoi[i_HI] = rhoi[i_HI]
                 - rhoi[i_OH]    * 1.0/17.0
                 - rhoi[i_H2O]   * 2.0/18.0
                 - rhoi[i_CH]    * 1.0/13.0
                 - rhoi[i_CH2]   * 2.0/14.0;
                  

  /* charge neutrality */
  rhoi[i_elec   ] = 0.0
          + rhoi[i_HII    ]
          + rhoi[i_HeII   ] / 4.0
          + rhoi[i_HeIII  ] / 4.0 * 2.0;
  if(grackle_data->primordial_chemistry > 1)
    rhoi[i_elec] = rhoi[i_elec]
          - rhoi[i_HM     ]
          + rhoi[i_H2II   ] / 2.0;
  if(grackle_data->primordial_chemistry > 2)
    rhoi[i_elec] = rhoi[i_elec]
          + rhoi[i_DII    ] / 2.0;
  if(grackle_data->primordial_chemistry > 3)
    rhoi[i_elec] = rhoi[i_elec]
          - rhoi[i_DM     ] / 2.0
          + rhoi[i_HDII   ] / 3.0
          + rhoi[i_HeHII  ] / 5.0;
  if(grackle_data->metal_chemistry)
    if(Metallicity > 1.0e-9) 
      rhoi[i_elec] = rhoi[i_elec]
            + rhoi[i_CII    ] / 12.0
            + rhoi[i_COII   ] / 28.0
            + rhoi[i_OII    ] / 16.0
            + rhoi[i_OHII   ] / 17.0
            + rhoi[i_H2OII  ] / 18.0
            + rhoi[i_H3OII  ] / 19.0
            + rhoi[i_O2II   ] / 32.0;

  *rho_gas = n0 * MPROTON / grackle_data->HydrogenFractionByMass;
//printf("DENSITY %13.5e  %13.5e  %13.5e \n", n0, MPROTON, grackle_data->HydrogenFractionByMass);

  nHmH = n0 * MPROTON;

  rhoi[i_elec   ]  *= *rho_gas * 0.00054462;
  rhoi[i_HI     ]  *= *rho_gas;
  rhoi[i_HII    ]  *= *rho_gas;
  rhoi[i_HeI    ]  *= *rho_gas;
  rhoi[i_HeII   ]  *= *rho_gas;
  rhoi[i_HeIII  ]  *= *rho_gas;
  rhoi[i_HM     ]  *= *rho_gas;
  rhoi[i_H2I    ]  *= *rho_gas;
  rhoi[i_H2II   ]  *= *rho_gas;
  rhoi[i_DI     ]  *= *rho_gas;
  rhoi[i_DII    ]  *= *rho_gas;
  rhoi[i_HDI    ]  *= *rho_gas;
  rhoi[i_HeHII  ]  *= *rho_gas;
  rhoi[i_DM     ]  *= *rho_gas;
  rhoi[i_HDII   ]  *= *rho_gas;
  rhoi[i_CI     ]  *= *rho_gas;
  rhoi[i_CII    ]  *= *rho_gas;
  rhoi[i_CO     ]  *= *rho_gas;
  rhoi[i_CO2    ]  *= *rho_gas;
  rhoi[i_OI     ]  *= *rho_gas;
  rhoi[i_OH     ]  *= *rho_gas;
  rhoi[i_H2O    ]  *= *rho_gas;
  rhoi[i_O2     ]  *= *rho_gas;
  rhoi[i_SiI    ]  *= *rho_gas;
  rhoi[i_SiOI   ]  *= *rho_gas;
  rhoi[i_SiO2I  ]  *= *rho_gas;
  rhoi[i_CH     ]  *= *rho_gas;
  rhoi[i_CH2    ]  *= *rho_gas;
  rhoi[i_COII   ]  *= *rho_gas;
  rhoi[i_OII    ]  *= *rho_gas;
  rhoi[i_OHII   ]  *= *rho_gas;
  rhoi[i_H2OII  ]  *= *rho_gas;
  rhoi[i_H3OII  ]  *= *rho_gas;
  rhoi[i_O2II   ]  *= *rho_gas;
  rhoi[i_Mg     ]  *= *rho_gas;
  rhoi[i_Al     ]  *= *rho_gas;
  rhoi[i_S      ]  *= *rho_gas;
  rhoi[i_Fe     ]  *= *rho_gas;
  rhoi[i_SiM    ]  *= *rho_gas;
  rhoi[i_FeM    ]  *= *rho_gas;
  rhoi[i_Mg2SiO4]  *= *rho_gas;
  rhoi[i_MgSiO3 ]  *= *rho_gas;
  rhoi[i_Fe3O4  ]  *= *rho_gas;
  rhoi[i_AC     ]  *= *rho_gas;
  rhoi[i_SiO2D  ]  *= *rho_gas;
  rhoi[i_MgO    ]  *= *rho_gas;
  rhoi[i_FeS    ]  *= *rho_gas;
  rhoi[i_Al2O3  ]  *= *rho_gas;
  rhoi[i_reforg  ] *= *rho_gas;
  rhoi[i_volorg  ] *= *rho_gas;
  rhoi[i_H2Oice  ] *= *rho_gas;

  /* floor */
//for(isp = 0; isp < 52; isp++)
//  if(rhoi[isp] / grackle_units.density_units < 1.0e-20)
//    rhoi[isp] = 1.0e-20 * grackle_units.density_units;

//printf("Z      %13.7f\n", Metallicity);
//printf("A(C)   %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XC  / 12.0 / grackle_data->HydrogenFractionByMass));
//printf("[C /H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XC  / 12.0 / grackle_data->HydrogenFractionByMass) - 8.43);
//printf("[O /H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XO  / 16.0 / grackle_data->HydrogenFractionByMass) - 8.69);
//printf("[Mg/H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XMg / 24.0 / grackle_data->HydrogenFractionByMass) - 7.60);
//printf("[Al/H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XAl / 27.0 / grackle_data->HydrogenFractionByMass) - 6.45);
//printf("[Si/H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XSi / 28.0 / grackle_data->HydrogenFractionByMass) - 7.51);
//printf("[S /H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XS  / 32.0 / grackle_data->HydrogenFractionByMass) - 7.12);
//printf("[Fe/H] %13.7f\n", 12.0 + log10(Metallicity * grackle_data->loc_XFe / 56.0 / grackle_data->HydrogenFractionByMass) - 7.50);

//printf("d  ");                              printf("%10.2e ", *rho_gas  / grackle_units.density_units); printf("\n");
//printf("p1 "); for(isp = 0; isp < 6; isp++) printf("%10.2e ", rhoi[i_SiM    ] / nHmH / 28.0); printf("\n");
//printf("p2 "); for(isp = 6; isp < 9; isp++) printf("%10.2e ", rhoi[i_FeM    ] / nHmH / 56.0); printf("\n");
//printf("p3 "); for(isp = 9; isp <12; isp++) printf("%10.2e ", rhoi[i_Mg2SiO4] / nHmH / 70.0); printf("\n");
//printf("p4 "); for(isp =12; isp <15; isp++) printf("%10.2e ", rhoi[i_MgSiO3 ] / nHmH /100.0); printf("\n");
//printf("C  "); for(isp =15; isp <19; isp++) printf("%10.2e ", rhoi[i_Fe3O4  ] / nHmH / 77.3); printf("\n");
//printf("O  "); for(isp =19; isp <23; isp++) printf("%10.2e ", rhoi[i_AC     ] / nHmH / 12.0); printf("\n");
//printf("Si "); for(isp =23; isp <26; isp++) printf("%10.2e ", rhoi[i_SiO2D  ] / nHmH / 60.0); printf("\n");
//printf("c  "); for(isp =26; isp <30; isp++) printf("%10.2e ", rhoi[i_MgO    ] / nHmH / 40.0); printf("\n");
//printf("o  "); for(isp =30; isp <34; isp++) printf("%10.2e ", rhoi[i_FeS    ] / nHmH / 88.0); printf("\n");
//printf("M  "); for(isp =34; isp <38; isp++) printf("%10.2e ", rhoi[i_Al2O3  ] / nHmH / 51.0); printf("\n");
//printf("d1 "); for(isp =38; isp <43; isp++) printf("%10.2e ", ); printf("\n");
//printf("d2 "); for(isp =43; isp <48; isp++) printf("%10.2e ", ); printf("\n");
//printf("d3 "); for(isp =48; isp <51; isp++) printf("%10.2e ", ); printf("\n");

  printf("d  ");
  fprintf(stdout, "%10.2e ", n0);
  fprintf(stdout, "%10.2e ", T0);
  printf("\n");
  printf("p1 ");
  fprintf(stdout, "%10.2e ", rhoi[i_elec   ] / nHmH /  0.00054462); 
  fprintf(stdout, "%10.2e ", rhoi[i_HI     ] / nHmH /  1.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HII    ] / nHmH /  1.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HeI    ] / nHmH /  4.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HeII   ] / nHmH /  4.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HeIII  ] / nHmH /  4.0);       
  printf("\n");
  printf("p2 ");
  fprintf(stdout, "%10.2e ", rhoi[i_HM     ] / nHmH /  1.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_H2I    ] / nHmH /  2.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_H2II   ] / nHmH /  2.0);       
  printf("\n");
  printf("p3 ");
  fprintf(stdout, "%10.2e ", rhoi[i_DI     ] / nHmH /  2.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_DII    ] / nHmH /  2.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HDI    ] / nHmH /  3.0);       
  printf("\n");
  printf("p4 ");
  fprintf(stdout, "%10.2e ", rhoi[i_DM     ] / nHmH /  2.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HDII   ] / nHmH /  3.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_HeHII  ] / nHmH /  5.0);       
  printf("\n");
  printf("C  ");
  fprintf(stdout, "%10.2e ", rhoi[i_CI     ] / nHmH / 12.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_CII    ] / nHmH / 12.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_CO     ] / nHmH / 28.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_CO2    ] / nHmH / 44.0);       
  printf("\n");
  printf("O  ");
  fprintf(stdout, "%10.2e ", rhoi[i_OI     ] / nHmH / 16.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_OH     ] / nHmH / 17.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_H2O    ] / nHmH / 18.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_O2     ] / nHmH / 32.0);       
  printf("\n");
  printf("Si ");
  fprintf(stdout, "%10.2e ", rhoi[i_SiI    ] / nHmH / 28.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_SiOI   ] / nHmH / 44.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_SiO2I  ] / nHmH / 60.0);       
  printf("\n");
  printf("Z1 ");
  fprintf(stdout, "%10.2e ", rhoi[i_CH     ] / nHmH / 13.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_CH2    ] / nHmH / 14.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_COII   ] / nHmH / 28.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_OII    ] / nHmH / 16.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_OHII   ] / nHmH / 17.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_H2OII  ] / nHmH / 18.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_H3OII  ] / nHmH / 19.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_O2II   ] / nHmH / 32.0);       
  printf("\n");
  printf("Z2 ");
  fprintf(stdout, "%10.2e ", rhoi[i_Mg     ] / nHmH / 24.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_Al     ] / nHmH / 27.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_S      ] / nHmH / 32.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_Fe     ] / nHmH / 56.0);       
  printf("\n");
  printf("d1 ");
  fprintf(stdout, "%10.2e ", rhoi[i_MgSiO3 ] / nHmH /100.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_AC     ] / nHmH / 12.0);       
  printf("\n");
  printf("d2 ");
  fprintf(stdout, "%10.2e ", rhoi[i_SiM    ] / nHmH / 28.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_FeM    ] / nHmH / 56.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_Mg2SiO4] / nHmH / 70.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_Fe3O4  ] / nHmH / 77.3);       
  fprintf(stdout, "%10.2e ", rhoi[i_SiO2D  ] / nHmH / 60.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_MgO    ] / nHmH / 40.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_FeS    ] / nHmH / 88.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_Al2O3  ] / nHmH / 51.0);       
  printf("\n");
  printf("d3 ");
  fprintf(stdout, "%10.2e ", rhoi[i_volorg ] / nHmH / 32.0);       
  fprintf(stdout, "%10.2e ", rhoi[i_reforg ] / nHmH / 22.68);       
  fprintf(stdout, "%10.2e ", rhoi[i_H2Oice ] / nHmH / 18.0);       
  printf("\n");

  *gamma   = get_ai(T0, rhoi);
  *energy  = *rho_gas/get_mu(rhoi) * KBOLTZ*T0 / (*gamma-1.0);
   printf("ENERGY %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e \n"
     , *rho_gas, get_mu(rhoi), KBOLTZ, T0, *gamma, *energy);

}


void update     (double dt, double *density, double *pressure, double *energy
               , double *edot_adia, double edot_cool, double *gamma, double *rhoi
               , double *specific_heat_ratio, double *dt_col, int itchem)
{
  int itr, isp;
  int n_sp = 6 + grackle_data->primordial_chemistry * 3;
  n_sp = 51;
  double t_col,t_col_a,t_col0;
  double Gamma,mu,rho,rho0,drho,e0,de,e,p0,dp,p,T;
  double G_comp,L_net,SHR;
  double err,err0;

  if(itchem==0) {
    return;
  }

  rho0   = *density;
  p0     = *pressure;
  e0     = *energy/rho0;
  G_comp = *edot_adia/rho0;
  L_net  =  edot_cool/rho0;
  Gamma  = *gamma;
  mu     =  get_mu(rhoi);
  SHR    = *specific_heat_ratio;
  t_col0 =  sqrt(3.0 * M_PI / (32.0 * GRAVITY * rho0));

  rho = rho0;
  t_col = t_col0;
  e = e0;

    drho = rho0/t_col*dt;
    G_comp = (Gamma-1.0) * e0/t_col;
#ifdef ISOTHERMAL
    de = 0.0;
#else
    de = (G_comp+L_net)*dt;
#endif
    rho = rho0+drho;
    e = e0+de;
    T = e*(Gamma-1.0)*mu/KBOLTZ;
    Gamma = get_ai(T,rhoi);
    p = rho*KBOLTZ*T/mu;
    SHR = (p-p0)/p*t_col/dt;

  *density = rho;
  *pressure = p;
  *energy = e*rho;
  *edot_adia = G_comp*rho;
  *gamma = Gamma;
  *specific_heat_ratio = SHR;
  *dt_col = t_col;

  for(isp=0; isp<n_sp; isp++) rhoi[isp] *= rho / rho0;

}



void update_phys(double dt, double *density, double *pressure, double *energy
               , double *edot_adia, double edot_cool, double *gamma, double *rhoi
               , double *specific_heat_ratio, double *dt_col, int itchem)
{
  int itr, isp;
  int n_sp = 6 + grackle_data->primordial_chemistry * 3;
  n_sp = 51;
  double t_col,t_col_a,t_col0;
  double Gamma,mu,rho,rho0,drho,e0,de,e,p0,dp,p,T;
  double G_comp,L_net,SHR;
  double err,err0;

//if(itchem==0) {
//  return;
//}

  rho0   = *density;
  p0     = *pressure;
  e0     = *energy/rho0;
  G_comp = *edot_adia/rho0;
  L_net  =  edot_cool/rho0;
  Gamma  = *gamma;
  mu     =  get_mu(rhoi);
  SHR    = *specific_heat_ratio;
  t_col0 =  compute_collapse_time(SHR,rho0);

  rho = rho0;
  t_col = t_col0;
  e = e0;
#ifndef GRAIN_GROWTH_TEST
  itr=0;
  err = 100.0;
  while(itr<=1000 && fabs(err)>1.0e-4) {
    drho = rho0/t_col*dt;
    G_comp = (Gamma-1.0) * e0/t_col;
#ifdef ISOTHERMAL
    de = 0.0;
#else
    de = (G_comp+L_net)*dt;
#endif
    rho = rho0+drho;
    e = e0+de;
    T = e*(Gamma-1.0)*mu/KBOLTZ;
    Gamma = get_ai(T,rhoi);
    p = rho*KBOLTZ*T/mu;
    SHR = (p-p0)/p*t_col/dt;
    t_col_a = compute_collapse_time(SHR,rho0);

    err = (t_col_a-t_col) / t_col_a;
    if(itr>0) t_col_a = t_col - err*(t_col-t_col0)/(err-err0);
    err0 = err;
    t_col0 = t_col;
    t_col  = t_col_a;
    itr++;
  }
#endif

  *density = rho;
  *pressure = p;
  *energy = e*rho;
  *edot_adia = G_comp*rho;
  *gamma = Gamma;
  *specific_heat_ratio = SHR;
  *dt_col = t_col;

  for(isp=0; isp<n_sp; isp++) rhoi[isp] *= rho / rho0;

}



double compute_collapse_time(double SHR, double rho_gas)
{
   double t_col0, t_col;
   double factor;

  t_col0 = sqrt(1.0/24.0/M_PI/GRAVITY/rho_gas);

#ifdef ISOTHERMAL
      factor = 0.6;
#else
  if(rho_gas < 1.0e2*MPROTON/grackle_data->HydrogenFractionByMass) {
      factor = 0.0;
  } else {
    if(SHR<0.83)
      factor = 0.0;
    else if(SHR<1.0)
      factor = 0.6 + 2.5*(SHR-1.0) - 6.0*pow(SHR-1.0,2.0);
    else if(SHR<4.0/3.0)
      factor = 1.0 + 0.2*(SHR-4.0/3.0) - 2.9*pow(SHR-4.0/3.0,2.0);
    else
      factor = 1.0;

    if(factor>0.95) factor = 0.95;
  }
#endif

  t_col = t_col0 / sqrt(1.0-factor);
  t_col *= 3.0*M_PI/2.0;
  t_col *= pow(10.0, -1.0 + 0.1*if0);
//      / pow(rho_gas / 1e-24, 1.0/16.0);

  return t_col;
}


double timestep(int *itchem, double dt_col, double dt_chem, double dt_cool, double dt)
{
  double dt_new;
  double Courant = 0.02;

  dt_chem = fabs(dt_chem);
  dt_cool = fabs(dt_cool);
//printf("[timestep] dt %9.2e dt_col %9.2e  dt_cool %9.2e  dt_chem %9.2e ",dt,dt_col,dt_cool,dt_chem);
  if(! *itchem) {
    dt_new = dt_chem;
    if(dt<dt_chem) *itchem = 1;
  } else {
    if(dt_col<dt_cool)
      dt_new = Courant*dt_col;
    else
      dt_new = Courant*dt_cool;

    if(dt_new>2.0*dt)
      dt_new = 2.0*dt;
  }
// dt_new = Courant * dt_col;

//printf("%13.5e\n", dt_col);


//printf("==> dt_phys %9.2e (%1d) \n",dt_new, *itchem);

  return dt_new;

}


void init_grackle()
{
//printf("%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e \n"
//  , grackle_data->A_C  
//  , grackle_data->A_O  
//  , grackle_data->A_Mg 
//  , grackle_data->A_Al 
//  , grackle_data->A_Si 
//  , grackle_data->A_S  
//  , grackle_data->A_Fe 
//);

  my_fields.grid_rank = 3;
  my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_start     = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_end       = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_dx = 0.0;
  my_fields.density          = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.internal_energy  = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.x_velocity       = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.y_velocity       = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.z_velocity       = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.       e_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      HI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     HII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     HeI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    HeII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   HeIII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      HM_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     H2I_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    H2II_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      DI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     DII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     HDI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      DM_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    HDII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   HeHII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      CI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     CII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      CO_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     CO2_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      OI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      OH_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     H2O_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      O2_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     SiI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    SiOI_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   SiO2I_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      CH_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     CH2_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    COII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     OII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    OHII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   H2OII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   H3OII_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.    O2II_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      Mg_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      Al_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.       S_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      Fe_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     SiM_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     FeM_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields. Mg2SiO4_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.  MgSiO3_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   Fe3O4_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.      AC_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   SiO2D_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     MgO_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.     FeS_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   Al2O3_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.  metal_density  = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.  reforg_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.  volorg_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.  H2Oice_density = malloc( ncell * ncell * ncell *sizeof(double));
  if (grackle_data->use_dust_density_field)
    my_fields.  dust_density = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_loc     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_C13     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_C20     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_C25     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_C30     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_F13     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_F15     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_F50     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_F80     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_P170    = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_P200    = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.   metal_Y19     = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.volumetric_heating_rate = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.specific_heating_rate   = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_heating_rate         = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.isrf_habing             = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_HI_ionization_rate   = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_HeI_ionization_rate  = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_HeII_ionization_rate = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_H2_dissociation_rate = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_HDI_dissociation_rate = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_CI_ionization_rate    = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_OI_ionization_rate    = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_CO_dissociation_rate  = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_OH_dissociation_rate  = malloc( ncell * ncell * ncell *sizeof(double));
  my_fields.RT_H2O_dissociation_rate = malloc( ncell * ncell * ncell *sizeof(double));
}


double get_mu(double rhoi[])
{

  return MPROTON * (
    rhoi[i_elec   ] 
  + rhoi[i_HI     ] 
  + rhoi[i_HII    ] 
  + rhoi[i_HeI    ] 
  + rhoi[i_HeII   ] 
  + rhoi[i_HeIII  ] 
  + rhoi[i_HM     ] 
  + rhoi[i_H2I    ] 
  + rhoi[i_H2II   ]
 ) / (
    rhoi[i_elec   ] / 0.00054462
  + rhoi[i_HI     ] / 1.0
  + rhoi[i_HII    ] / 1.0
  + rhoi[i_HeI    ] / 4.0
  + rhoi[i_HeII   ] / 4.0
  + rhoi[i_HeIII  ] / 4.0 
  + rhoi[i_HM     ] / 1.0
  + rhoi[i_H2I    ] / 2.0
  + rhoi[i_H2II   ] / 2.0 
  );

}

double get_ai(double temp, double rhoi[])
{
  double x, gamma_H2_minus1_inv;

  x=6100.0/temp;
  if(x>100) gamma_H2_minus1_inv = 2.5;
  else      gamma_H2_minus1_inv = 0.5*(5.0 + 2.0*x*x * exp(x)/(exp(x)-1.0)/(exp(x)-1.0));

  return 1.0 + 
  ( rhoi[i_elec   ] / 0.00054462
  + rhoi[i_HI     ] / 1.0
  + rhoi[i_HII    ] / 1.0
  + rhoi[i_HeI    ] / 4.0
  + rhoi[i_HeII   ] / 4.0
  + rhoi[i_HeIII  ] / 4.0
  + rhoi[i_HM     ] / 1.0
  + rhoi[i_H2I    ] / 2.0
  + rhoi[i_H2II   ] / 2.0
 ) / (
                    1.5 * rhoi[i_elec   ] / 0.00054462
  +                 1.5 * rhoi[i_HI     ] / 1.0
  +                 1.5 * rhoi[i_HII    ] / 1.0
  +                 1.5 * rhoi[i_HeI    ] / 4.0
  +                 1.5 * rhoi[i_HeII   ] / 4.0
  +                 1.5 * rhoi[i_HeIII  ] / 4.0
  +                 1.5 * rhoi[i_HM     ] / 1.0
  + gamma_H2_minus1_inv * rhoi[i_H2I    ] / 2.0
  + gamma_H2_minus1_inv * rhoi[i_H2II   ] / 2.0
 );

}

int compute_abundances(int itime, double dt, double rho, double *e, double *rhoi)
{
  int ip, isp;
//printf("before\n");
//printf("p1 "); for(isp = 0; isp < 6; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p2 "); for(isp = 6; isp < 9; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p3 "); for(isp = 9; isp <12; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p4 "); for(isp =12; isp <15; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("C  "); for(isp =15; isp <19; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("O  "); for(isp =19; isp <23; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("Si "); for(isp =23; isp <26; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("c  "); for(isp =26; isp <30; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("o  "); for(isp =30; isp <34; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("M  "); for(isp =34; isp <38; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("d1 "); for(isp =38; isp <43; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("d2 "); for(isp =43; isp <48; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("d3 "); for(isp =48; isp <51; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("e  "); printf("%10.2e\n", *e);

  double dtFixed = dt / grackle_units.time_units;
  my_fields.grid_dimension[0] = ncell;
  my_fields.grid_dimension[1] = ncell;
  my_fields.grid_dimension[2] = ncell;
  my_fields.grid_start[0] = 0; my_fields.grid_end[0] = ncell - 1;
  my_fields.grid_start[1] = 0; my_fields.grid_end[1] = ncell - 1;
  my_fields.grid_start[2] = 0; my_fields.grid_end[2] = ncell - 1;

  for(ip = 0; ip < ncell * ncell * ncell ; ip++) {
  my_fields.density          [ip] = rho / grackle_units.density_units;
  my_fields.internal_energy  [ip] = *e/rho / (pow(grackle_units.a_units * grackle_units.length_units / grackle_units.time_units, 2));
  my_fields.x_velocity       [ip] = 0.0;
  my_fields.y_velocity       [ip] = 0.0;
  my_fields.z_velocity       [ip] = 0.0;
  my_fields.       e_density [ip] = rhoi[i_elec   ] / 0.00054462 / grackle_units.density_units;
  my_fields.      HI_density [ip] = rhoi[i_HI     ] / grackle_units.density_units;
  my_fields.     HII_density [ip] = rhoi[i_HII    ] / grackle_units.density_units;
  my_fields.     HeI_density [ip] = rhoi[i_HeI    ] / grackle_units.density_units;
  my_fields.    HeII_density [ip] = rhoi[i_HeII   ] / grackle_units.density_units;
  my_fields.   HeIII_density [ip] = rhoi[i_HeIII  ] / grackle_units.density_units;
  my_fields.      HM_density [ip] = rhoi[i_HM     ] / grackle_units.density_units;
  my_fields.     H2I_density [ip] = rhoi[i_H2I    ] / grackle_units.density_units;
  my_fields.    H2II_density [ip] = rhoi[i_H2II   ] / grackle_units.density_units;
  my_fields.      DI_density [ip] = rhoi[i_DI     ] / grackle_units.density_units;
  my_fields.     DII_density [ip] = rhoi[i_DII    ] / grackle_units.density_units;
  my_fields.     HDI_density [ip] = rhoi[i_HDI    ] / grackle_units.density_units;
  my_fields.      DM_density [ip] = rhoi[i_DM     ] / grackle_units.density_units;
  my_fields.    HDII_density [ip] = rhoi[i_HDII   ] / grackle_units.density_units;
  my_fields.   HeHII_density [ip] = rhoi[i_HeHII  ] / grackle_units.density_units;
  my_fields.      CI_density [ip] = rhoi[i_CI     ] / grackle_units.density_units;
  my_fields.     CII_density [ip] = rhoi[i_CII    ] / grackle_units.density_units;
  my_fields.      CO_density [ip] = rhoi[i_CO     ] / grackle_units.density_units;
  my_fields.     CO2_density [ip] = rhoi[i_CO2    ] / grackle_units.density_units;
  my_fields.      OI_density [ip] = rhoi[i_OI     ] / grackle_units.density_units;
  my_fields.      OH_density [ip] = rhoi[i_OH     ] / grackle_units.density_units;
  my_fields.     H2O_density [ip] = rhoi[i_H2O    ] / grackle_units.density_units;
  my_fields.      O2_density [ip] = rhoi[i_O2     ] / grackle_units.density_units;
  my_fields.     SiI_density [ip] = rhoi[i_SiI    ] / grackle_units.density_units;
  my_fields.    SiOI_density [ip] = rhoi[i_SiOI   ] / grackle_units.density_units;
  my_fields.   SiO2I_density [ip] = rhoi[i_SiO2I  ] / grackle_units.density_units;
  my_fields.      CH_density [ip] = rhoi[i_CH     ] / grackle_units.density_units;
  my_fields.     CH2_density [ip] = rhoi[i_CH2    ] / grackle_units.density_units;
  my_fields.    COII_density [ip] = rhoi[i_COII   ] / grackle_units.density_units;
  my_fields.     OII_density [ip] = rhoi[i_OII    ] / grackle_units.density_units;
  my_fields.    OHII_density [ip] = rhoi[i_OHII   ] / grackle_units.density_units;
  my_fields.   H2OII_density [ip] = rhoi[i_H2OII  ] / grackle_units.density_units;
  my_fields.   H3OII_density [ip] = rhoi[i_H3OII  ] / grackle_units.density_units;
  my_fields.    O2II_density [ip] = rhoi[i_O2II   ] / grackle_units.density_units;
  my_fields.      Mg_density [ip] = rhoi[i_Mg     ] / grackle_units.density_units;
  my_fields.      Al_density [ip] = rhoi[i_Al     ] / grackle_units.density_units;
  my_fields.       S_density [ip] = rhoi[i_S      ] / grackle_units.density_units;
  my_fields.      Fe_density [ip] = rhoi[i_Fe     ] / grackle_units.density_units;
  my_fields.     SiM_density [ip] = rhoi[i_SiM    ] / grackle_units.density_units;
  my_fields.     FeM_density [ip] = rhoi[i_FeM    ] / grackle_units.density_units;
  my_fields. Mg2SiO4_density [ip] = rhoi[i_Mg2SiO4] / grackle_units.density_units;
  my_fields.  MgSiO3_density [ip] = rhoi[i_MgSiO3 ] / grackle_units.density_units;
  my_fields.   Fe3O4_density [ip] = rhoi[i_Fe3O4  ] / grackle_units.density_units;
  my_fields.      AC_density [ip] = rhoi[i_AC     ] / grackle_units.density_units;
  my_fields.   SiO2D_density [ip] = rhoi[i_SiO2D  ] / grackle_units.density_units;
  my_fields.     MgO_density [ip] = rhoi[i_MgO    ] / grackle_units.density_units;
  my_fields.     FeS_density [ip] = rhoi[i_FeS    ] / grackle_units.density_units;
  my_fields.   Al2O3_density [ip] = rhoi[i_Al2O3  ] / grackle_units.density_units;
  my_fields.  reforg_density [ip] = rhoi[i_reforg  ] / grackle_units.density_units;
  my_fields.  volorg_density [ip] = rhoi[i_volorg  ] / grackle_units.density_units;
  my_fields.  H2Oice_density [ip] = rhoi[i_H2Oice  ] / grackle_units.density_units;
  my_fields.    metal_density[ip] = Metallicity * my_fields.density[ip];
  if(grackle_data->multi_metals) {
  my_fields.    metal_loc    [ip] = Metallicity_loc * my_fields.density[ip];
  my_fields.    metal_C13    [ip] = Metallicity_C13 * my_fields.density[ip];
  my_fields.    metal_C20    [ip] = Metallicity_C20 * my_fields.density[ip];
  my_fields.    metal_C25    [ip] = Metallicity_C25 * my_fields.density[ip];
  my_fields.    metal_C30    [ip] = Metallicity_C30 * my_fields.density[ip];
  my_fields.    metal_F13    [ip] = Metallicity_F13 * my_fields.density[ip];
  my_fields.    metal_F15    [ip] = Metallicity_F15 * my_fields.density[ip];
  my_fields.    metal_F50    [ip] = Metallicity_F50 * my_fields.density[ip];
  my_fields.    metal_F80    [ip] = Metallicity_F80 * my_fields.density[ip];
  my_fields.    metal_P170   [ip] = Metallicity_P170* my_fields.density[ip];
  my_fields.    metal_P200   [ip] = Metallicity_P200* my_fields.density[ip];
  my_fields.    metal_Y19    [ip] = Metallicity_Y19 * my_fields.density[ip];
  my_fields.    metal_density[ip] = my_fields.metal_loc [ip]
                                  + my_fields.metal_C13 [ip]
                                  + my_fields.metal_C20 [ip]
                                  + my_fields.metal_C25 [ip]
                                  + my_fields.metal_C30 [ip]
                                  + my_fields.metal_F13 [ip]
                                  + my_fields.metal_F15 [ip]
                                  + my_fields.metal_F50 [ip]
                                  + my_fields.metal_F80 [ip]
                                  + my_fields.metal_P170[ip]
                                  + my_fields.metal_P200[ip]
                                  + my_fields.metal_Y19 [ip]; 
  }
  if (grackle_data->use_dust_density_field)
    my_fields.  dust_density[ip] = 0.0;

  double l_jeans, NH, Av, G0;
  double kdiss_H2I, kdiss_HDI, kdiss_CO;
  if(test == 2) {
    l_jeans = sqrt(M_PI * (5.0/3.0) * (5.0/3.0 - 1.0) * *e/rho / (GRAVITY * rho));
    NH = (XH * rho / MPROTON) * l_jeans;
//  printf("%13.5e %13.5e\n",(XH * rho / MPROTON), l_jeans);
    if(Metallicity < 1.0e-8 * Zsun)
      Av = 0.0;
    else
      Av = 5.3e-22 * NH * (Metallicity / Zsun);
//  printf("%13.7f\n", Av);
    if(iD0 < 9)
      G0 = pow(10.0, -iD0) * G00;
    else
      G0 = 0.0;
 
    kdiss_H2I = 4.50e-11 * exp(-2.5 * Av) * G0;
    kdiss_HDI = 4.50e-11 * exp(-2.5 * Av) * G0;
    kdiss_CO  = 2.0e-10  * exp(-3.5 * Av) * G0 / 1.7;
  }

//if(grackle_data->radiative_transfer_use_H2_shielding) {
//  kdiss_H2I *= rhoi[i_H2I] / (2.0 * MPROTON) * l_jeans
//}

//printf("%13.5e %13.5e %13.5e\n", (XH * rho / MPROTON), (l_jeans / 3.0856e18), G);

//double N_H2I, N_HDI, N_CO;
//double x, b_doppler;
//N_H2I = rhoi[i_H2I] / ( 2.0 * MPROTON) * l_jeans;
//x = N_H2I / 5.0e14;
//b_doppler = 1E-5_DKIND *
//        sqrt(2._DKIND * kboltz *
//             tgas1d(i) / (2._DKIND * mass_h))
//f_shield = 0.965_DKIND / (1._DKIND + x/b_doppler)**1.1+
//  0.035_DKIND * exp(-8.5E-4_DKIND * sqrt(1._DKIND +x))/
//  sqrt(1._DKIND + x)

//N_HDI = rhoi[i_HDI] / ( 3.0 * MPROTON) * l_jeans;
//N_CO  = rhoi[i_CO ] / (28.0 * MPROTON) * l_jeans;

  my_fields.volumetric_heating_rate[ip]  = 0.0;
  my_fields.specific_heating_rate  [ip]  = 0.0;
  my_fields.RT_heating_rate        [ip]  = 0.0;
  my_fields.RT_HI_ionization_rate  [ip]  = 0.0;
  my_fields.RT_HeI_ionization_rate [ip]  = 0.0;
  my_fields.RT_HeII_ionization_rate[ip]  = 0.0;
  my_fields.RT_H2_dissociation_rate[ip]  = kdiss_H2I * grackle_units.time_units;
  my_fields.RT_HDI_dissociation_rate[ip] = kdiss_HDI * grackle_units.time_units;
  my_fields.RT_CI_ionization_rate   [ip] = 0.0;
  my_fields.RT_OI_ionization_rate   [ip] = 0.0;
  my_fields.RT_CO_dissociation_rate [ip] = kdiss_CO  * grackle_units.time_units;
  my_fields.RT_OH_dissociation_rate [ip] = 0.0;
  my_fields.RT_H2O_dissociation_rate[ip] = 0.0;
//printf("# %13.5e %13.5e\n", my_fields.RT_H2_dissociation_rate[ip], grackle_units.time_units);

  /* Compute BH radiation field */
  double Mbh = 100.0; // Msun
  double alpha = 0.5;
  double nu0 = 10.0;    // eV
  double nu1 = 100.0e3; // eV
  double l, n0, T0;
  double Ledd = 3.3e6 * LSOLAR * (Mbh / 100.0);
  double L = 1.e-2 * Ledd;
  double F, radius[ncell];

  if(test == 4 && grackle_data->use_isrf_field) {
    /* Flux */
//  radius[ip] = pow(10.0, -5.0 + 0.1 * ip) * CM_PER_PC;
    radius[ip] = 3e-4 * CM_PER_PC;
    F = L / (4.0 * M_PI * pow(radius[ip],2)); // erg / s / cm^2
    my_fields.isrf_habing[ip] = F / 5.3e-3; // cgs -> habing unit
    if(itime == 0) {
      printf("Flux %13.5e cm %13.5e erg / s / cm^2\n", radius[ip], F);
      printf("Flux %13.5e Habing\n", my_fields.isrf_habing[ip]);
    }
    my_fields.RT_HI_ionization_rate   [ip] = F / E_HI   * sigma_HI   * grackle_units.time_units;
    my_fields.RT_HeI_ionization_rate  [ip] = F / E_HeI  * sigma_HeI  * grackle_units.time_units;
    my_fields.RT_HeII_ionization_rate [ip] = F / E_HeII * sigma_HeII * grackle_units.time_units;
    my_fields.RT_H2_dissociation_rate [ip] = F / E_H2I  * sigma_H2I  * grackle_units.time_units;
    my_fields.RT_heating_rate         [ip] = F          * sigma_HI;  // cgs
  }

  }
  
//#define TEST_OUTPUT
#ifdef TEST_OUTPUT
  for(ip = 0; ip < ncell * ncell * ncell ; ip++) {
    printf("        d %13.5e\n", my_fields.density        [ip]);
    printf("        e %13.5e\n", my_fields.internal_energy[ip]);
    printf("       de %13.5e\n", my_fields.      e_density[ip]);
    printf("       HI %13.5e\n", my_fields.     HI_density[ip]);
    printf("      HII %13.5e\n", my_fields.    HII_density[ip]);
    printf("      HeI %13.5e\n", my_fields.    HeI_density[ip]);
    printf("     HeII %13.5e\n", my_fields.   HeII_density[ip]);
    printf("    HeIII %13.5e\n", my_fields.  HeIII_density[ip]);
    printf("       HM %13.5e\n", my_fields.     HM_density[ip]);
    printf("      H2I %13.5e\n", my_fields.    H2I_density[ip]);
    printf("     H2II %13.5e\n", my_fields.   H2II_density[ip]);
    printf("       DI %13.5e\n", my_fields.     DI_density[ip]);
    printf("      DII %13.5e\n", my_fields.    DII_density[ip]);
    printf("      HDI %13.5e\n", my_fields.    HDI_density[ip]);
    printf("       DM %13.5e\n", my_fields.     DM_density[ip]);
    printf("     HDII %13.5e\n", my_fields.   HDII_density[ip]);
    printf("    HeHII %13.5e\n", my_fields.  HeHII_density[ip]);
    printf("       CI %13.5e\n", my_fields.     CI_density[ip]);
    printf("      CII %13.5e\n", my_fields.    CII_density[ip]);
    printf("       CO %13.5e\n", my_fields.     CO_density[ip]);
    printf("      CO2 %13.5e\n", my_fields.    CO2_density[ip]);
    printf("       OI %13.5e\n", my_fields.     OI_density[ip]);
    printf("       OH %13.5e\n", my_fields.     OH_density[ip]);
    printf("      H2O %13.5e\n", my_fields.    H2O_density[ip]);
    printf("       O2 %13.5e\n", my_fields.     O2_density[ip]);
    printf("      SiI %13.5e\n", my_fields.    SiI_density[ip]);
    printf("     SiOI %13.5e\n", my_fields.   SiOI_density[ip]);
    printf("    SiO2I %13.5e\n", my_fields.  SiO2I_density[ip]);
    printf("       CH %13.5e\n", my_fields.     CH_density[ip]);
    printf("      CH2 %13.5e\n", my_fields.    CH2_density[ip]);
    printf("     COII %13.5e\n", my_fields.   COII_density[ip]);
    printf("      OII %13.5e\n", my_fields.    OII_density[ip]);
    printf("     OHII %13.5e\n", my_fields.   OHII_density[ip]);
    printf("    H2OII %13.5e\n", my_fields.  H2OII_density[ip]);
    printf("    H3OII %13.5e\n", my_fields.  H3OII_density[ip]);
    printf("     O2II %13.5e\n", my_fields.   O2II_density[ip]);
    printf("       Mg %13.5e\n", my_fields.     Mg_density[ip]);
    printf("       Al %13.5e\n", my_fields.     Al_density[ip]);
    printf("        S %13.5e\n", my_fields.      S_density[ip]);
    printf("       Fe %13.5e\n", my_fields.     Fe_density[ip]);
    printf("      SiM %13.5e\n", my_fields.    SiM_density[ip]);
    printf("      FeM %13.5e\n", my_fields.    FeM_density[ip]);
    printf("  Mg2SiO4 %13.5e\n", my_fields.Mg2SiO4_density[ip]);
    printf("   MgSiO3 %13.5e\n", my_fields. MgSiO3_density[ip]);
    printf("    Fe3O4 %13.5e\n", my_fields.  Fe3O4_density[ip]);
    printf("       AC %13.5e\n", my_fields.     AC_density[ip]);
    printf("    SiO2D %13.5e\n", my_fields.  SiO2D_density[ip]);
    printf("      MgO %13.5e\n", my_fields.    MgO_density[ip]);
    printf("      FeS %13.5e\n", my_fields.    FeS_density[ip]);
    printf("    Al2O3 %13.5e\n", my_fields.  Al2O3_density[ip]);
    printf("   reforg %13.5e\n", my_fields.  reforg_density[ip]);
    printf("   volorg %13.5e\n", my_fields.  volorg_density[ip]);
    printf("   H2Oice %13.5e\n", my_fields.  H2Oice_density[ip]);
    printf("    metal %13.5e\n", my_fields.  metal_density[ip] / my_fields.density        [ip]);
    printf("     dust %13.5e\n", my_fields.   dust_density[ip]);
    printf("metal loc  %13.5e\n", my_fields.metal_loc [ip]);
    printf("metal C13  %13.5e\n", my_fields.metal_C13 [ip]);
    printf("metal C20  %13.5e\n", my_fields.metal_C20 [ip]);
    printf("metal C25  %13.5e\n", my_fields.metal_C25 [ip]);
    printf("metal C30  %13.5e\n", my_fields.metal_C30 [ip]);
    printf("metal F13  %13.5e\n", my_fields.metal_F13 [ip]);
    printf("metal F15  %13.5e\n", my_fields.metal_F15 [ip]);
    printf("metal F50  %13.5e\n", my_fields.metal_F50 [ip]);
    printf("metal F80  %13.5e\n", my_fields.metal_F80 [ip]);
    printf("metal P170 %13.5e\n", my_fields.metal_P170[ip]);
    printf("metal P200 %13.5e\n", my_fields.metal_P200[ip]);
    printf("metal Y19  %13.5e\n", my_fields.metal_Y19 [ip]);
  }
#endif


//printf("\nSolve chemistry\n");
  if (solve_chemistry(&grackle_units,
                      &my_fields,
                      (double)dtFixed
                    ) == FAIL) {
    printf("Error in Grackle solve_chemistry.\n");
    return FAIL;
  }

  rhoi[i_elec   ] = my_fields.       e_density [0] * 0.00054462 * grackle_units.density_units;
  rhoi[i_HI     ] = my_fields.      HI_density [0] * grackle_units.density_units;
  rhoi[i_HII    ] = my_fields.     HII_density [0] * grackle_units.density_units;
  rhoi[i_HeI    ] = my_fields.     HeI_density [0] * grackle_units.density_units;
  rhoi[i_HeII   ] = my_fields.    HeII_density [0] * grackle_units.density_units;
  rhoi[i_HeIII  ] = my_fields.   HeIII_density [0] * grackle_units.density_units;
  rhoi[i_HM     ] = my_fields.      HM_density [0] * grackle_units.density_units;
  rhoi[i_H2I    ] = my_fields.     H2I_density [0] * grackle_units.density_units;
  rhoi[i_H2II   ] = my_fields.    H2II_density [0] * grackle_units.density_units;
  rhoi[i_DI     ] = my_fields.      DI_density [0] * grackle_units.density_units;
  rhoi[i_DII    ] = my_fields.     DII_density [0] * grackle_units.density_units;
  rhoi[i_HDI    ] = my_fields.     HDI_density [0] * grackle_units.density_units;
  rhoi[i_DM     ] = my_fields.      DM_density [0] * grackle_units.density_units;
  rhoi[i_HDII   ] = my_fields.    HDII_density [0] * grackle_units.density_units;
  rhoi[i_HeHII  ] = my_fields.   HeHII_density [0] * grackle_units.density_units;
  rhoi[i_CI     ] = my_fields.      CI_density [0] * grackle_units.density_units;
  rhoi[i_CII    ] = my_fields.     CII_density [0] * grackle_units.density_units;
  rhoi[i_CO     ] = my_fields.      CO_density [0] * grackle_units.density_units;
  rhoi[i_CO2    ] = my_fields.     CO2_density [0] * grackle_units.density_units;
  rhoi[i_OI     ] = my_fields.      OI_density [0] * grackle_units.density_units;
  rhoi[i_OH     ] = my_fields.      OH_density [0] * grackle_units.density_units;
  rhoi[i_H2O    ] = my_fields.     H2O_density [0] * grackle_units.density_units;
  rhoi[i_O2     ] = my_fields.      O2_density [0] * grackle_units.density_units;
  rhoi[i_SiI    ] = my_fields.     SiI_density [0] * grackle_units.density_units;
  rhoi[i_SiOI   ] = my_fields.    SiOI_density [0] * grackle_units.density_units;
  rhoi[i_SiO2I  ] = my_fields.   SiO2I_density [0] * grackle_units.density_units;
  rhoi[i_CH     ] = my_fields.      CH_density [0] * grackle_units.density_units;
  rhoi[i_CH2    ] = my_fields.     CH2_density [0] * grackle_units.density_units;
  rhoi[i_COII   ] = my_fields.    COII_density [0] * grackle_units.density_units;
  rhoi[i_OII    ] = my_fields.     OII_density [0] * grackle_units.density_units;
  rhoi[i_OHII   ] = my_fields.    OHII_density [0] * grackle_units.density_units;
  rhoi[i_H2OII  ] = my_fields.   H2OII_density [0] * grackle_units.density_units;
  rhoi[i_H3OII  ] = my_fields.   H3OII_density [0] * grackle_units.density_units;
  rhoi[i_O2II   ] = my_fields.    O2II_density [0] * grackle_units.density_units;
  rhoi[i_Mg     ] = my_fields.      Mg_density [0] * grackle_units.density_units;
  rhoi[i_Al     ] = my_fields.      Al_density [0] * grackle_units.density_units;
  rhoi[i_S      ] = my_fields.       S_density [0] * grackle_units.density_units;
  rhoi[i_Fe     ] = my_fields.      Fe_density [0] * grackle_units.density_units;
  rhoi[i_SiM    ] = my_fields.     SiM_density [0] * grackle_units.density_units;
  rhoi[i_FeM    ] = my_fields.     FeM_density [0] * grackle_units.density_units;
  rhoi[i_Mg2SiO4] = my_fields. Mg2SiO4_density [0] * grackle_units.density_units;
  rhoi[i_MgSiO3 ] = my_fields.  MgSiO3_density [0] * grackle_units.density_units;
  rhoi[i_Fe3O4  ] = my_fields.   Fe3O4_density [0] * grackle_units.density_units;
  rhoi[i_AC     ] = my_fields.      AC_density [0] * grackle_units.density_units;
  rhoi[i_SiO2D  ] = my_fields.   SiO2D_density [0] * grackle_units.density_units;
  rhoi[i_MgO    ] = my_fields.     MgO_density [0] * grackle_units.density_units;
  rhoi[i_FeS    ] = my_fields.     FeS_density [0] * grackle_units.density_units;
  rhoi[i_Al2O3  ] = my_fields.   Al2O3_density [0] * grackle_units.density_units;
  rhoi[i_reforg  ] = my_fields.  reforg_density [0] * grackle_units.density_units;
  rhoi[i_volorg  ] = my_fields.  volorg_density [0] * grackle_units.density_units;
  rhoi[i_H2Oice  ] = my_fields.  H2Oice_density [0] * grackle_units.density_units;
  *e = my_fields.internal_energy[0]*rho * (pow(grackle_units.a_units * grackle_units.length_units / grackle_units.time_units, 2));

//printf("DD %13.5e %13.5e \n", rhoi[i_HI]/rho, rhoi[i_OH]/rho);
//printf("after\n");
//printf("p1 "); for(isp = 0; isp < 6; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p2 "); for(isp = 6; isp < 9; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p3 "); for(isp = 9; isp <12; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("p4 "); for(isp =12; isp <15; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("C  "); for(isp =15; isp <19; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("O  "); for(isp =19; isp <23; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("Si "); for(isp =23; isp <26; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("c  "); for(isp =26; isp <30; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("o  "); for(isp =30; isp <34; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("M  "); for(isp =34; isp <38; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("d1 "); for(isp =38; isp <43; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("d2 "); for(isp =43; isp <48; isp++) printf("%10.2e ", rhoi[isp] / rho); printf("\n");
//printf("e  "); printf("%10.2e\n", *e);
  return SUCCESS;

}


void output(int itime, double time, double rho, double energy, double rhoi[], double edot_adia)
{
  int iter, isp;
  int n_sp = 6 + grackle_data->primordial_chemistry * 3;
  n_sp = 51;
  double nH, T, y_H2I;
  double gamma0, gamma, red;
  double mu;
  double nHmH;
  double charge, H_tot, D_tot, He_tot;
  double Z_tot, C_tot, O_tot, Mg_tot, Al_tot, Si_tot, S_tot, Fe_tot;
  double H_frac, D_frac, He_frac;
  double C_frac, O_frac, Mg_frac, Al_frac, Si_frac, S_frac, Fe_frac;
  double SiM_frac     , FeM_frac    , Mg2SiO4_frac, MgSiO3_frac , Fe3O4_frac  
       , AC_frac      , SiO2D_frac  , MgO_frac    , FeS_frac    , Al2O3_frac
       , reforg_frac  , volorg_frac , H2Oice_frac ;

  nH = grackle_data->HydrogenFractionByMass * rho / MPROTON;

  y_H2I = rhoi[i_H2I] / rho / grackle_data->HydrogenFractionByMass / 2.0;
  if(y_H2I<0.25) gamma0 = 5.0/3.0;
  else           gamma0 = 1.4;

  iter=0;
  do {
    T = (gamma0-1.0) * energy / rho * get_mu(rhoi) / KBOLTZ;
    gamma = get_ai(T, rhoi);
    red = fabs(gamma - gamma0)/gamma0;
    gamma0 = gamma;
    iter++;
  } while (red>1.0e-10 && iter<100);
//printf("%13.5e %13.5e %13.5e %13.5e\n", nH, get_mu(rhoi)/MPROTON, T, gamma);

  nHmH = grackle_data->HydrogenFractionByMass * rho;

   H_frac = (1.0 - grackle_data->DeuteriumToHydrogenRatio) * grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
// H_frac = grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
   D_frac = grackle_data->DeuteriumToHydrogenRatio * grackle_data->HydrogenFractionByMass * (1.0 - Metallicity);
  He_frac = (1.0 - grackle_data->HydrogenFractionByMass) * (1.0 - Metallicity);

  if(Metallicity > 1.0e-9) {
//  if(grackle_data->grain_growth) {
          C_frac = grackle_data->SN0_fC [grackle_data->metal_abundances] * Metallicity;
          O_frac = grackle_data->SN0_fO [grackle_data->metal_abundances] * Metallicity;
         Mg_frac = grackle_data->SN0_fMg[grackle_data->metal_abundances] * Metallicity;
         Al_frac = grackle_data->SN0_fAl[grackle_data->metal_abundances] * Metallicity;
         Si_frac = grackle_data->SN0_fSi[grackle_data->metal_abundances] * Metallicity;
          S_frac = grackle_data->SN0_fS [grackle_data->metal_abundances] * Metallicity;
         Fe_frac = grackle_data->SN0_fFe[grackle_data->metal_abundances] * Metallicity;
        SiM_frac = grackle_data->SN0_fSiM    [grackle_data->metal_abundances] * Metallicity;
        FeM_frac = grackle_data->SN0_fFeM    [grackle_data->metal_abundances] * Metallicity;
    Mg2SiO4_frac = grackle_data->SN0_fMg2SiO4[grackle_data->metal_abundances] * Metallicity;
     MgSiO3_frac = grackle_data->SN0_fMgSiO3 [grackle_data->metal_abundances] * Metallicity;
      Fe3O4_frac = grackle_data->SN0_fFe3O4  [grackle_data->metal_abundances] * Metallicity;
         AC_frac = grackle_data->SN0_fAC     [grackle_data->metal_abundances] * Metallicity;
      SiO2D_frac = grackle_data->SN0_fSiO2D  [grackle_data->metal_abundances] * Metallicity;
        MgO_frac = grackle_data->SN0_fMgO    [grackle_data->metal_abundances] * Metallicity;
        FeS_frac = grackle_data->SN0_fFeS    [grackle_data->metal_abundances] * Metallicity;
      Al2O3_frac = grackle_data->SN0_fAl2O3  [grackle_data->metal_abundances] * Metallicity;
     reforg_frac = grackle_data->SN0_freforg [grackle_data->metal_abundances] * Metallicity;
     volorg_frac = grackle_data->SN0_fvolorg [grackle_data->metal_abundances] * Metallicity;
     H2Oice_frac = grackle_data->SN0_fH2Oice [grackle_data->metal_abundances] * Metallicity;
//  } else {
//        C_frac = grackle_data->SN0_XC [grackle_data->metal_abundances] * Metallicity;
//        O_frac = grackle_data->SN0_XO [grackle_data->metal_abundances] * Metallicity;
//       Mg_frac = grackle_data->SN0_XMg[grackle_data->metal_abundances] * Metallicity;
//       Al_frac = grackle_data->SN0_XAl[grackle_data->metal_abundances] * Metallicity;
//       Si_frac = grackle_data->SN0_XSi[grackle_data->metal_abundances] * Metallicity;
//        S_frac = grackle_data->SN0_XS [grackle_data->metal_abundances] * Metallicity;
//       Fe_frac = grackle_data->SN0_XFe[grackle_data->metal_abundances] * Metallicity;
//  }
  } else {
   C_frac = 0.0;
   O_frac = 0.0;
  Mg_frac = 0.0;
  Al_frac = 0.0;
  Si_frac = 0.0;
   S_frac = 0.0;
  Fe_frac = 0.0;
  }

    charge =
          - rhoi[i_elec   ] / 0.00054462
          + rhoi[i_HII    ]
          + rhoi[i_HeII   ] / 4.0
          + rhoi[i_HeIII  ] / 4.0 * 2.0;
  if(grackle_data->primordial_chemistry > 1)
    charge = charge
          - rhoi[i_HM     ]
          + rhoi[i_H2II   ] / 2.0;
  if(grackle_data->primordial_chemistry > 2)
    charge = charge
          + rhoi[i_DII    ] / 2.0;
  if(grackle_data->primordial_chemistry > 3)
    charge = charge
          - rhoi[i_DM     ] / 2.0
          + rhoi[i_HDII   ] / 3.0
          + rhoi[i_HeHII  ] / 5.0;
  if(grackle_data->metal_chemistry)
    if(Metallicity > 1.0e-9)
      charge = charge
            + rhoi[i_CII    ] / 12.0
            + rhoi[i_COII   ] / 28.0
            + rhoi[i_OII    ] / 16.0
            + rhoi[i_OHII   ] / 17.0
            + rhoi[i_H2OII  ] / 18.0
            + rhoi[i_H3OII  ] / 19.0
            + rhoi[i_O2II   ] / 32.0;
//charge = charge / (H_frac * rho);


    H_tot = 
            rhoi[i_HI     ]
          + rhoi[i_HII    ];
  if(grackle_data->primordial_chemistry > 1)
    H_tot = H_tot
          + rhoi[i_HM     ]
          + rhoi[i_H2I    ] 
          + rhoi[i_H2II   ];
  if(grackle_data->primordial_chemistry > 2)
    H_tot = H_tot
          + rhoi[i_HDI    ] /  3.0;
  if(grackle_data->primordial_chemistry > 3)
    H_tot = H_tot
          + rhoi[i_HDII   ] /  3.0
          + rhoi[i_HeHII  ] /  5.0;
  if(grackle_data->metal_chemistry)
    if(Metallicity > 1.0e-9) 
      H_tot = H_tot
            + rhoi[i_OH     ] / 17.0
            + rhoi[i_H2O    ] / 18.0 * 2.0
            + rhoi[i_CH     ] / 13.0
            + rhoi[i_CH2    ] / 14.0 * 2.0
            + rhoi[i_OHII   ] / 17.0
            + rhoi[i_H2OII  ] / 18.0 * 2.0
            + rhoi[i_H3OII  ] / 19.0 * 3.0;
//H_tot = H_tot / (H_frac * rho) - 1.0;


    D_tot = 0.0;
  if(grackle_data->primordial_chemistry > 2)
    D_tot = D_tot
          + rhoi[i_DI     ]
          + rhoi[i_DII    ]
          + rhoi[i_HDI    ] /  3.0 * 2.0;
  if(grackle_data->primordial_chemistry > 3)
    D_tot = D_tot
          + rhoi[i_DM     ]
          + rhoi[i_HDII   ] /  3.0 * 2.0;
//D_tot = D_tot / (D_frac * rho) - 1.0;


    He_tot = 
            rhoi[i_HeI    ]
          + rhoi[i_HeII   ]
          + rhoi[i_HeIII  ];
  if(grackle_data->primordial_chemistry > 3)
    He_tot = He_tot
          + rhoi[i_HeHII  ] /  5.0 * 4.0;
//He_tot = He_tot / (He_frac * rho) - 1.0;


  if(Metallicity > 1.0e-9)
  if(grackle_data->metal_chemistry) {
 
//  Z_tot = my_fields.metal_density[0] * grackle_units.density_units;
    Z_tot = Metallicity * rho;

    C_tot = rhoi[i_CI     ]
          + rhoi[i_CII    ]
          + rhoi[i_CO     ] / 28.0 * 12.0
          + rhoi[i_CO2    ] / 44.0 * 12.0
          + rhoi[i_CH     ] / 13.0 * 12.0
          + rhoi[i_CH2    ] / 14.0 * 12.0
          + rhoi[i_COII   ] / 28.0 * 12.0;
    if(grackle_data->grain_growth) {
    C_tot = C_tot
          + rhoi[i_AC     ]
          + rhoi[i_reforg ] / 22.68 * 12.0
          + rhoi[i_volorg ] / 32.0  * 12.0;
    }
//C_tot = C_tot / (C_frac * rho) - 1.0;


     
    O_tot = rhoi[i_CO     ] / 28.0 * 16.0
          + rhoi[i_CO2    ] / 44.0 * 32.0
          + rhoi[i_OI     ]
          + rhoi[i_OH     ] / 17.0 * 16.0
          + rhoi[i_H2O    ] / 18.0 * 16.0
          + rhoi[i_O2     ]
          + rhoi[i_SiOI   ] / 44.0 * 16.0
          + rhoi[i_SiO2I  ] / 60.0 * 32.0
          + rhoi[i_COII   ] / 28.0 * 12.0
          + rhoi[i_OII    ]
          + rhoi[i_OHII   ] / 17.0 * 16.0
          + rhoi[i_H2OII  ] / 18.0 * 16.0
          + rhoi[i_H3OII  ] / 19.0 * 16.0
          + rhoi[i_O2II   ];
    if(grackle_data->grain_growth) {
    O_tot = O_tot
          + rhoi[i_Mg2SiO4] /140.0 * 64.0
          + rhoi[i_MgSiO3 ] /100.0 * 48.0
          + rhoi[i_Fe3O4  ] /232.0 * 64.0
          + rhoi[i_SiO2D  ] / 60.0 * 32.0
          + rhoi[i_MgO    ] / 40.0 * 16.0
          + rhoi[i_Al2O3  ] /102.0 * 48.0
          + rhoi[i_reforg ] / 22.68 * 8.0
          + rhoi[i_volorg ] / 32.0  * 16.0
          + rhoi[i_H2Oice ] / 18.0  * 16.0;
    }
//O_tot = O_tot / (O_frac * rho) - 1.0;


    Mg_tot = rhoi[i_Mg     ];
    if(grackle_data->grain_growth) {
    Mg_tot = Mg_tot
           + rhoi[i_Mg2SiO4] /140.0 * 48.0
           + rhoi[i_MgSiO3 ] /100.0 * 24.0
           + rhoi[i_MgO    ] / 40.0 * 24.0;
    }
//Mg_tot = Mg_tot / (Mg_frac * rho) - 1.0;


    Al_tot = rhoi[i_Al     ];
    if(grackle_data->grain_growth) {
    Al_tot = Al_tot
           + rhoi[i_Al2O3  ]/102.0 * 54.0;
    }
//Al_tot = Al_tot / (Al_frac * rho) - 1.0;


    Si_tot= rhoi[i_SiI    ]
          + rhoi[i_SiOI   ]/ 44.0 * 28.0
          + rhoi[i_SiO2I  ]/ 60.0 * 28.0;
    if(grackle_data->grain_growth) {
    Si_tot = Si_tot
          + rhoi[i_SiM    ]
          + rhoi[i_Mg2SiO4]/140.0 * 28.0
          + rhoi[i_MgSiO3 ]/100.0 * 28.0
          + rhoi[i_SiO2D  ]/ 60.0 * 28.0;
    }
//Si_tot = Si_tot / (Si_frac * rho) - 1.0;


    S_tot = rhoi[i_S      ];
    if(grackle_data->grain_growth) {
    S_tot = S_tot
          + rhoi[i_FeS    ]/ 88.0 * 32.0;
    }
//S_tot = S_tot / (S_frac * rho) - 1.0;


    Fe_tot = rhoi[i_Fe     ];
    if(grackle_data->grain_growth) {
    Fe_tot = Fe_tot
           + rhoi[i_FeM    ]
           + rhoi[i_Fe3O4  ]/232.0 *168.0
           + rhoi[i_FeS    ]/ 88.0 * 56.0;
    }
//Fe_tot = Fe_tot / (Fe_frac * rho) - 1.0;

  }

/*
  if(itime       == 0) {
    printf("%13s", "nH");
    printf("%13s", "charge");
    printf("%13s", " H_tot");
    printf("%13s", " D_tot");
    printf("%13s", "He_tot");
    printf("%13s", " Z_tot");
    printf("%13s", " C_tot");
    printf("%13s", " O_tot");
    printf("%13s", "Mg_tot");
    printf("%13s", "Al_tot");
    printf("%13s", "Si_tot");
    printf("%13s", " S_tot");
    printf("%13s", "Fe_tot");
    printf("\n");

    printf("%13.5e", 0.0);
    printf("%13.5e", 0.0);
    printf("%13.5e",  H_frac);
    printf("%13.5e",  D_frac);
    printf("%13.5e", He_frac);
    printf("%13.5e", Metallicity);
    printf("%13.5e",  C_frac);
    printf("%13.5e",  O_frac);
    printf("%13.5e", Mg_frac);
    printf("%13.5e", Al_frac);
    printf("%13.5e", Si_frac);
    printf("%13.5e",  S_frac);
    printf("%13.5e", Fe_frac);
    printf("\n");
  }
  if(itime % 100 == 0) {
    printf("%13.5e", nH);
    printf("%13.5e", charge / rho);
    printf("%13.5e",  H_tot / rho);
    printf("%13.5e",  D_tot / rho);
    printf("%13.5e", He_tot / rho);
    printf("%13.5e",  Z_tot / rho);
    printf("%13.5e",  C_tot / rho);
    printf("%13.5e",  O_tot / rho);
    printf("%13.5e", Mg_tot / rho);
    printf("%13.5e", Al_tot / rho);
    printf("%13.5e", Si_tot / rho);
    printf("%13.5e",  S_tot / rho);
    printf("%13.5e", Fe_tot / rho);
    printf("\n");
  }
*/

//printf("\nCalculate dust temperature\n");
  gr_float *dust_temperature    ;
  dust_temperature     = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.SiM_temperature      = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.FeM_temperature      = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.Mg2SiO4_temperature  = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.MgSiO3_temperature   = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.Fe3O4_temperature    = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.AC_temperature       = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.SiO2D_temperature    = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.MgO_temperature      = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.FeS_temperature      = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.Al2O3_temperature    = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.reforg_temperature   = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.volorg_temperature   = malloc(ncell * ncell * ncell * sizeof(double));
  my_fields.H2Oice_temperature   = malloc(ncell * ncell * ncell * sizeof(double));

  if (calculate_dust_temperature(&grackle_units,
                                 &my_fields
                                , dust_temperature    
                               ) == FAIL) {
    printf("Error in Grackle calculate_dust_temperature.\n");
    exit(1);
  }
//printf("TTT3 %13.5e\n", my_fields.MgSiO3_temperature[0]);

#ifdef GRAIN_GROWTH_TEST
//printf("%13.5e %13.5e\n", nH, T);
  fprintf(fdout, "%23.15e %13.5e ", time, T);
#else
  fprintf(fdout, "%13.5e %13.5e ", nH, T);
#endif
//fprintf(fdout, "%13.5e %13.5e ", nH, energy / rho);
  fprintf(fdout, "%13.5e ", rhoi[i_elec   ] / nHmH /  0.00054462); //  3 
  fprintf(fdout, "%13.5e ", rhoi[i_HI     ] / nHmH /  1.0);        //  4
  fprintf(fdout, "%13.5e ", rhoi[i_HII    ] / nHmH /  1.0);        //  5
  fprintf(fdout, "%13.5e ", rhoi[i_H2I    ] / nHmH /  2.0);        //  6
  fprintf(fdout, "%13.5e ", rhoi[i_HM     ] / nHmH /  1.0);        //  7
  fprintf(fdout, "%13.5e ", rhoi[i_H2II   ] / nHmH /  2.0);        //  8
  fprintf(fdout, "%13.5e ", rhoi[i_HeHII  ] / nHmH /  5.0);        //  9
  fprintf(fdout, "%13.5e ", rhoi[i_HeI    ] / nHmH /  4.0);        // 10
  fprintf(fdout, "%13.5e ", rhoi[i_HeII   ] / nHmH /  4.0);        // 11
  fprintf(fdout, "%13.5e ", rhoi[i_HeIII  ] / nHmH /  4.0);        // 12
  fprintf(fdout, "%13.5e ", rhoi[i_DI     ] / nHmH /  2.0);        // 13
  fprintf(fdout, "%13.5e ", rhoi[i_DII    ] / nHmH /  2.0);        // 14
  fprintf(fdout, "%13.5e ", rhoi[i_DM     ] / nHmH /  2.0);        // 15
  fprintf(fdout, "%13.5e ", rhoi[i_HDI    ] / nHmH /  3.0);        // 16
  fprintf(fdout, "%13.5e ", rhoi[i_HDII   ] / nHmH /  3.0);        // 17
  fprintf(fdout, "%13.5e ", rhoi[i_CI     ] / nHmH / 12.0);        // 18
  fprintf(fdout, "%13.5e ", rhoi[i_CII    ] / nHmH / 12.0);        // 19
  fprintf(fdout, "%13.5e ", rhoi[i_CO     ] / nHmH / 28.0);        // 20
  fprintf(fdout, "%13.5e ", rhoi[i_CO2    ] / nHmH / 44.0);        // 21
  fprintf(fdout, "%13.5e ", rhoi[i_OI     ] / nHmH / 16.0);        // 22
  fprintf(fdout, "%13.5e ", rhoi[i_OH     ] / nHmH / 17.0);        // 23
  fprintf(fdout, "%13.5e ", rhoi[i_H2O    ] / nHmH / 18.0);        // 24
  fprintf(fdout, "%13.5e ", rhoi[i_O2     ] / nHmH / 32.0);        // 25
  fprintf(fdout, "%13.5e ", rhoi[i_SiI    ] / nHmH / 28.0);        // 26
  fprintf(fdout, "%13.5e ", rhoi[i_SiOI   ] / nHmH / 44.0);        // 27
  fprintf(fdout, "%13.5e ", rhoi[i_SiO2I  ] / nHmH / 60.0);        // 28
  fprintf(fdout, "%13.5e ", rhoi[i_CH     ] / nHmH / 13.0);        // 29
  fprintf(fdout, "%13.5e ", rhoi[i_CH2    ] / nHmH / 14.0);        // 30
  fprintf(fdout, "%13.5e ", rhoi[i_COII   ] / nHmH / 28.0);        // 31
  fprintf(fdout, "%13.5e ", rhoi[i_OII    ] / nHmH / 16.0);        // 32
  fprintf(fdout, "%13.5e ", rhoi[i_OHII   ] / nHmH / 17.0);        // 33
  fprintf(fdout, "%13.5e ", rhoi[i_H2OII  ] / nHmH / 18.0);        // 34
  fprintf(fdout, "%13.5e ", rhoi[i_H3OII  ] / nHmH / 19.0);        // 35
  fprintf(fdout, "%13.5e ", rhoi[i_O2II   ] / nHmH / 32.0);        // 36
  fprintf(fdout, "%13.5e ", rhoi[i_Mg     ] / nHmH / 24.0);        // 37
  fprintf(fdout, "%13.5e ", rhoi[i_Al     ] / nHmH / 27.0);        // 38
  fprintf(fdout, "%13.5e ", rhoi[i_S      ] / nHmH / 32.0);        // 39
  fprintf(fdout, "%13.5e ", rhoi[i_Fe     ] / nHmH / 56.0);        // 40
  fprintf(fdout, "%13.5e ", rhoi[i_SiM    ] / nHmH / 28.0);        // 41
  fprintf(fdout, "%13.5e ", rhoi[i_FeM    ] / nHmH / 56.0);        // 42
  fprintf(fdout, "%13.5e ", rhoi[i_Mg2SiO4] / nHmH / 70.0);        // 43
  fprintf(fdout, "%13.5e ", rhoi[i_MgSiO3 ] / nHmH /100.0);        // 44
  fprintf(fdout, "%13.5e ", rhoi[i_Fe3O4  ] / nHmH / 77.3);        // 45
  fprintf(fdout, "%13.5e ", rhoi[i_AC     ] / nHmH / 12.0);        // 46
  fprintf(fdout, "%13.5e ", rhoi[i_SiO2D  ] / nHmH / 60.0);        // 47
  fprintf(fdout, "%13.5e ", rhoi[i_MgO    ] / nHmH / 40.0);        // 48
  fprintf(fdout, "%13.5e ", rhoi[i_FeS    ] / nHmH / 88.0);        // 49
  fprintf(fdout, "%13.5e ", rhoi[i_Al2O3  ] / nHmH / 51.0);        // 50
  fprintf(fdout, "%13.5e ", rhoi[i_reforg ] / nHmH / 22.68);       // 51
  fprintf(fdout, "%13.5e ", rhoi[i_volorg ] / nHmH / 32.0);        // 52
  fprintf(fdout, "%13.5e ", rhoi[i_H2Oice ] / nHmH / 18.0);        // 53
  fprintf(fdout, "%13.5e",           dust_temperature    [0]);     // 54
  fprintf(fdout, "%13.5e", my_fields.SiM_temperature     [0]);     // 55
  fprintf(fdout, "%13.5e", my_fields.FeM_temperature     [0]);     // 56
  fprintf(fdout, "%13.5e", my_fields.Mg2SiO4_temperature [0]);     // 57
  fprintf(fdout, "%13.5e", my_fields.MgSiO3_temperature  [0]);     // 58
  fprintf(fdout, "%13.5e", my_fields.Fe3O4_temperature   [0]);     // 69
  fprintf(fdout, "%13.5e", my_fields.AC_temperature      [0]);     // 60
  fprintf(fdout, "%13.5e", my_fields.SiO2D_temperature   [0]);     // 61
  fprintf(fdout, "%13.5e", my_fields.MgO_temperature     [0]);     // 62
  fprintf(fdout, "%13.5e", my_fields.FeS_temperature     [0]);     // 63
  fprintf(fdout, "%13.5e", my_fields.Al2O3_temperature   [0]);     // 64
  fprintf(fdout, "%13.5e", my_fields.reforg_temperature  [0]);     // 65
  fprintf(fdout, "%13.5e", my_fields.volorg_temperature  [0]);     // 66
  fprintf(fdout, "%13.5e", my_fields.H2Oice_temperature  [0]);     // 67
  fprintf(fdout, "%13.5e ",             edot_adia / rho);          // 68
//fprintf(fdout, "%13.5e ", charge);                               // 51
//fprintf(fdout, "%13.5e ", H_tot);                                // 52
//fprintf(fdout, "%13.5e ", D_tot);                                // 53
//fprintf(fdout, "%13.5e ", He_tot);                               // 54
//fprintf(fdout, "%13.5e ", C_tot);                                // 55
//fprintf(fdout, "%13.5e ", O_tot);                                // 56
//fprintf(fdout, "%13.5e ", Mg_tot);                               // 57
//fprintf(fdout, "%13.5e ", Al_tot);                               // 58
//fprintf(fdout, "%13.5e ", Si_tot);                               // 59
//fprintf(fdout, "%13.5e ", S_tot);                                // 60
//fprintf(fdout, "%13.5e ", Fe_tot);                               // 61
  fprintf(fdout, "\n");
  fflush(fdout);

  free(dust_temperature    );

}
