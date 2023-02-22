#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>
#include <omp.h>

using namespace std;

double d_edin[9] = { 1,0,0, 0,1,0, 0,0,1 };
double true_parameters[6] = {4.085, -2.960, 1.08, 1.32, 0.97, 0.51};
double true_params[2] = {0.22, -0.36};
double true_parameter = 0.497;
double Array[6] = {};
double Array2[12] = {};
double E_cohesion;
int ind_min_global;
double a0_global;
double h;
double Params_bounds[6][2] = {{0.0, 0.1}, {0.0685, 0.1370}, {0.7853, 1.570666}, {7.2853, 14.5706}, {2.0927, 4.1853}, {1.9257, 3.8514}};

class Atom
{
    public:
        string name;
        double x, y, z;

        Atom(double _x, double _y, double _z, string _name)
        {
            setAtom(_x, _y, _z, _name);
        }

        void setAtom(double _x, double _y, double _z, string _name)
        {
            x = _x;
            y = _y;
            z = _z;
            name = _name;
        }

        void getAtom()
        {
            cout << name << endl;
            cout << x << endl;
            cout << y << endl;
            cout << z << endl;
            cout << endl;
        }

        void changeAtom(string _name)
        {
            name = _name;
        }
};

vector<Atom> GCK(string name, double a0)
{
    vector <Atom>Atoms_GCK;

    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Atoms_GCK.push_back(Atom((0+i)*a0, (0+j)*a0, (0+k)*a0, name));
                Atoms_GCK.push_back(Atom((0+i)*a0, (0.5+j)*a0, (0.5+k)*a0, name));
                Atoms_GCK.push_back(Atom((0.5+i)*a0, (0.5+j)*a0, (0+k)*a0, name));
                Atoms_GCK.push_back(Atom((0.5+i)*a0, (0+j)*a0, (0.5+k)*a0, name));
            }
        }
    }

    return Atoms_GCK;
}

double Energy(const vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0)
{
    double Energy = 0;

    for (int i=0; i<Vect.size(); i++)
    {
        double Er = 0, Eb = 0;

        double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
        double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
        double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

        double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

        for (int j=0; j<Vect.size(); j++)
        {
            for (int dx=-1; dx<2; dx++)
                for (int dy=-1; dy<2; dy++)
                    for (int dz=-1; dz<2; dz++)
                    {
                        if (i != j || dx !=0 || dy !=0 || dz !=0 )
                        {
                            double tmp_x = Vect[j].x + a0*3*dx;
                            double tmp_y = Vect[j].y + a0*3*dy;
                            double tmp_z = Vect[j].z + a0*3*dz;

                            tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                            tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                            tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                            tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                            tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                            tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                            double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                            if (rij <= cutoff)
                            {
                                Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                                Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                            }
                        }
                    }
        }

        Eb = -sqrt(Eb);
        Energy += Er + Eb;
    }

    return Energy;
}

double Energy_parallel(const vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, int i)
{
    double Energy = 0;

    double Er = 0, Eb = 0;

    double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
    double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
    double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

    double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

    for (int j=0; j<Vect.size(); j++)
    {
        for (int dx=-1; dx<2; dx++)
            for (int dy=-1; dy<2; dy++)
                for (int dz=-1; dz<2; dz++)
                {
                    if (i != j || dx !=0 || dy !=0 || dz !=0 )
                    {
                        double tmp_x = Vect[j].x + a0*3*dx;
                        double tmp_y = Vect[j].y + a0*3*dy;
                        double tmp_z = Vect[j].z + a0*3*dz;

                        tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                        tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                        tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                        tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                        tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                        tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                        double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                        if (rij <= cutoff)
                        {
                            Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                            Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                        }
                    }
                }
    }

    Eb = -sqrt(Eb);
    Energy += Er + Eb;

    return Energy;
}

double Energy_AB(vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[6])
{
    double Energy = 0;

    for (int i=0; i<Vect.size(); i++)
    {
        double Er = 0, Eb = 0;

        double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
        double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
        double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

        double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

        for (int j=0; j<Vect.size(); j++)
        {
            for (int dx=-1; dx<2; dx++)
                for (int dy=-1; dy<2; dy++)
                    for (int dz=-1; dz<2; dz++)
                    {
                        if (i != j || dx !=0 || dy !=0 || dz !=0 )
                        {
                            double tmp_x = Vect[j].x + a0*3*dx;
                            double tmp_y = Vect[j].y + a0*3*dy;
                            double tmp_z = Vect[j].z + a0*3*dz;

                            tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                            tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                            tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                            tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                            tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                            tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                            double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                            if (rij <= cutoff)
                            {
                                if((Vect[i].name == "Ag") && (Vect[j].name == "Ag"))
                                {
                                    Er = Er + ((x[0]/x[5])*(rij-x[5])+x[1])*exp(-x[3]*(rij/x[5]-1));

                                    Eb = Eb + x[2]*x[2]*exp(-2*x[4]*(rij/x[5]-1));
                                }
                                else if(((Vect[i].name == "V") && (Vect[j].name == "Ag")) || ((Vect[i].name == "Ag") && (Vect[j].name == "V")))
                                {
                                    Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                                    Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                                }
                            }
                        }
                    }
        }

        Eb = -sqrt(Eb);
        Energy += Er + Eb;
    }

    return Energy;
}

double Energy_AB_parallel(vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[6], int i)
{
    double Energy = 0;

    double rij_min = 10 * a0;

    double Er = 0, Eb = 0;

    double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
    double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
    double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

    double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

    for (int j=0; j<Vect.size(); j++)
    {
        for (int dx=-1; dx<2; dx++)
            for (int dy=-1; dy<2; dy++)
                for (int dz=-1; dz<2; dz++)
                {
                    if (i != j || dx !=0 || dy !=0 || dz !=0 )
                    {
                        double tmp_x = Vect[j].x + a0*3*dx;
                        double tmp_y = Vect[j].y + a0*3*dy;
                        double tmp_z = Vect[j].z + a0*3*dz;

                        tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                        tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                        tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                        tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                        tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                        tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                        double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                        if (rij <= cutoff)
                        {
                            if((Vect[i].name == "Ag") && (Vect[j].name == "Ag"))
                            {
                                Er = Er + ((x[0]/x[5])*(rij-x[5])+x[1])*exp(-x[3]*(rij/x[5]-1));

                                Eb = Eb + x[2]*x[2]*exp(-2*x[4]*(rij/x[5]-1));
                            }
                            else if(((Vect[i].name == "V") && (Vect[j].name == "Ag")) || ((Vect[i].name == "Ag") && (Vect[j].name == "V")))
                            {
                                Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                                Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                            }
                        }
                    }
                }
    }

    Eb = -sqrt(Eb);
    Energy += Er + Eb;

    return Energy;
}

double Energy_last(vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[12])
{
    double Energy = 0;

    for (int i=0; i<Vect.size(); i++)
    {
        double Er = 0, Eb = 0;

        double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
        double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
        double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

        double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

        for (int j=0; j<Vect.size(); j++)
        {
            for (int dx=-1; dx<2; dx++)
                for (int dy=-1; dy<2; dy++)
                    for (int dz=0; dz<1; dz++)
                    {
                        if (i != j || dx !=0 || dy !=0 || dz !=0 )
                        {
                            double tmp_x = Vect[j].x + a0*3*dx;
                            double tmp_y = Vect[j].y + a0*3*dy;
                            double tmp_z = Vect[j].z + a0*3*dz;

                            tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                            tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                            tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                            tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                            tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                            tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                            double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                            if (rij <= cutoff)
                            {
                                if((Vect[i].name == "Ag") && (Vect[j].name == "Ag"))
                                {
                                    Er = Er + ((x[0]/x[5])*(rij-x[5])+x[1])*exp(-x[3]*(rij/x[5]-1));

                                    Eb = Eb + x[2]*x[2]*exp(-2*x[4]*(rij/x[5]-1));
                                }
                                else if(((Vect[i].name == "V") && (Vect[j].name == "Ag")) || ((Vect[i].name == "Ag") && (Vect[j].name == "V")))
                                {
                                    Er = Er + ((x[6]/x[11])*(rij-x[11])+x[7])*exp(-x[9]*(rij/x[11]-1));

                                    Eb = Eb + x[8]*x[8]*exp(-2*x[10]*(rij/x[11]-1));
                                }
                                else if((Vect[i].name == "V") && (Vect[j].name == "V"))
                                {
                                    Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                                    Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                                }
                            }
                        }
                    }
        }

        Eb = -sqrt(Eb);
        Energy += Er + Eb;
    }

    return Energy;
}

double Energy_last_parallel(vector <Atom>&Vect, double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[12], int i)
{
    double Energy = 0;

    double Er = 0, Eb = 0;

    double a0_x = a0 * d_trans[0] + a0 * d_trans[1] + a0 * d_trans[2];
    double a0_y = a0 * d_trans[3] + a0 * d_trans[4] + a0 * d_trans[5];
    double a0_z = a0 * d_trans[6] + a0 * d_trans[7] + a0 * d_trans[8];

    double cutoff = 1.7 * fmax(fmax(a0_x, a0_y), a0_z);

    for (int j=0; j<Vect.size(); j++)
    {
        for (int dx=-1; dx<2; dx++)
            for (int dy=-1; dy<2; dy++)
                for (int dz=0; dz<1; dz++)
                {
                    if (i != j || dx !=0 || dy !=0 || dz !=0 )
                    {
                        double tmp_x = Vect[j].x + a0*3*dx;
                        double tmp_y = Vect[j].y + a0*3*dy;
                        double tmp_z = Vect[j].z + a0*3*dz;

                        tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                        tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                        tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];

                        tmp_x -= Vect[i].x * d_trans[0] + Vect[i].y * d_trans[1] + Vect[i].z * d_trans[2];
                        tmp_y -= Vect[i].x * d_trans[3] + Vect[i].y * d_trans[4] + Vect[i].z * d_trans[5];
                        tmp_z -= Vect[i].x * d_trans[6] + Vect[i].y * d_trans[7] + Vect[i].z * d_trans[8];

                        double rij = sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z);

                        if (rij <= cutoff)
                        {
                            if((Vect[i].name == "Ag") && (Vect[j].name == "Ag"))
                            {
                                Er = Er + ((x[0]/x[5])*(rij-x[5])+x[1])*exp(-x[3]*(rij/x[5]-1));

                                Eb = Eb + x[2]*x[2]*exp(-2*x[4]*(rij/x[5]-1));
                            }
                            else if(((Vect[i].name == "V") && (Vect[j].name == "Ag")) || ((Vect[i].name == "Ag") && (Vect[j].name == "V")))
                            {
                                Er = Er + ((x[6]/x[11])*(rij-x[11])+x[7])*exp(-x[9]*(rij/x[11]-1));

                                Eb = Eb + x[8]*x[8]*exp(-2*x[10]*(rij/x[11]-1));
                            }
                            else if((Vect[i].name == "V") && (Vect[j].name == "V"))
                            {
                                Er = Er + ((A1/r0)*(rij-r0)+A0)*exp(-p*(rij/r0-1));

                                Eb = Eb + ksi*ksi*exp(-2*q*(rij/r0-1));
                            }
                        }
                    }
                }
    }

    Eb = -sqrt(Eb);
    Energy += Er + Eb;

    return Energy;
}

double EnergyFinal(double rij, double x[6])
{
    double Energy, Er, Eb;

    Er = ((x[0]/x[5])*(rij-x[5])+x[1])*exp(-x[3]*(rij/x[5]-1));

    Eb = x[2]*x[2]*exp(-2*x[4]*(rij/x[5]-1));

    Eb = -sqrt(Eb);

    Energy = Er + Eb;

    return Energy;
}

double EnergySol(double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[6], double E_coh)
{
    vector <Atom> Vect = GCK("Ag", a0);

    Vect[0].changeAtom("V");

    //double E_AB = Energy_AB(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    double Energy_per_thread[4] = {0.0, 0.0, 0.0, 0.0};
    omp_set_dynamic(0);

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num] += Energy_AB_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    double E_AB = 0;

    for(int i=0; i<4; i++)
    {
        E_AB += Energy_per_thread[i];
    }

    double E_coh_A = -5.31;

    double E_coh_B = E_coh;

    double E_B = E_coh_B * Vect.size();

    Vect[0].changeAtom("Ag");

    return E_AB - E_B - E_coh_A + E_coh_B;
}

double EnergyIn(double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[12])
{
    int N = 2;

    vector <Atom> Vect = GCK("Ag", a0);

    /*Vect[81].changeAtom("V");
    Vect[83].changeAtom("V");

    double E_dim_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect[81].changeAtom("Ag");
    Vect[83].changeAtom("Ag");

    double E_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect[81].changeAtom("V");

    double E_adatom_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect[81].changeAtom("Ag");*/

    double Energy_per_thread[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    omp_set_dynamic(0);

    Vect[81].changeAtom("V");
    Vect[83].changeAtom("V");

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][0] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect[81].changeAtom("Ag");
    Vect[83].changeAtom("Ag");

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][1] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect[81].changeAtom("V");

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][2] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect[81].changeAtom("Ag");

    double E_dim_surf = 0, E_surf = 0, E_adatom_surf = 0;

    for(int i=0; i<4; i++)
    {
        E_dim_surf += Energy_per_thread[i][0];
        E_surf += Energy_per_thread[i][1];
        E_adatom_surf += Energy_per_thread[i][2];
    }

    return E_dim_surf - E_surf - N * (E_adatom_surf - E_surf);
}

double EnergyOn(double a0, double d_trans[9], double A1, double A0, double ksi, double p, double q, double r0, double x[12])
{
    int N = 2;

    vector <Atom> Vect = GCK("Ag", a0);

    /*Vect.push_back(Atom(3*a0-1.5*a0, 0.5*a0, 3*a0, "V"));
    Vect.push_back(Atom(2*a0, 0.0, 3*a0, "V"));

    double E_dim_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect.pop_back();
    Vect.pop_back();

    double E_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect.push_back(Atom(3*a0-1.5*a0, 0.5*a0, 3*a0, "V"));

    double E_adatom_surf = Energy_last(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x);

    Vect.pop_back();*/

    double Energy_per_thread[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    omp_set_dynamic(0);

    Vect.push_back(Atom(3*a0-1.5*a0, 0.5*a0, 3*a0, "V"));
    Vect.push_back(Atom(2*a0, 0.0, 3*a0, "V"));

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][0] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect.pop_back();
    Vect.pop_back();

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][1] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect.push_back(Atom(3*a0-1.5*a0, 0.5*a0, 3*a0, "V"));

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][2] += Energy_last_parallel(Vect, a0, d_trans, A1, A0, ksi, p, q, r0, x, i);
        }
    }

    Vect.pop_back();

    double E_dim_surf = 0, E_surf = 0, E_adatom_surf = 0;

    for(int i=0; i<4; i++)
    {
        E_dim_surf += Energy_per_thread[i][0];
        E_surf += Energy_per_thread[i][1];
        E_adatom_surf += Energy_per_thread[i][2];
    }

    return E_dim_surf - E_surf - N * (E_adatom_surf - E_surf);
}

void Parameters(double &E, double a0, double &B, double &c11, double &c12, double &c44, double A1, double A0, double ksi, double p, double q, double r0)
{
    double v0 = a0*a0*a0/4, const_p = 0.8018993929636421, alpha = 0.001, alpha2 = 0.000001;

    double d_b_plus[9]={1+alpha, 0,0,0, 1+alpha, 0,0,0, 1+alpha};

    double d_b_minus[9]={1-alpha, 0,0,0, 1-alpha, 0,0,0, 1-alpha};

    double d_c11_plus[9]={1+alpha, 0,0,0, 1+alpha, 0,0,0, 1};

    double d_c11_minus[9]={1-alpha, 0,0,0, 1-alpha, 0,0,0, 1};

    double d_c12_plus[9]={1+alpha, 0,0,0, 1-alpha, 0,0,0, 1};

    double d_c12_minus[9]={1-alpha, 0,0,0, 1+alpha, 0,0,0, 1};

    double d_c44_plus[9]={1, alpha, 0, alpha, 1, 0, 0,0, 1/(1-alpha2)};

    double d_c44_minus[9]={1, -alpha, 0, -alpha, 1, 0, 0,0, 1/(1-alpha2)};

    vector <Atom> Vect_p = GCK("Ag", a0);

    vector <Atom> Vect_B_plus = GCK("Ag", a0);

    vector <Atom> Vect_B_minus = GCK("Ag", a0);

    vector <Atom> Vect_c11_plus = GCK("Ag", a0);

    vector <Atom> Vect_c11_minus = GCK("Ag", a0);

    vector <Atom> Vect_c12_plus = GCK("Ag", a0);

    vector <Atom> Vect_c12_minus = GCK("Ag", a0);

    vector <Atom> Vect_c44_plus = GCK("Ag", a0);

    vector <Atom> Vect_c44_minus = GCK("Ag", a0);

    /*double E_p = Energy(Vect_p, a0, d_edin, A1, A0, ksi, p, q, r0) / Vect_p.size();

    double E_B_plus = Energy(Vect_B_plus, a0, d_b_plus, A1, A0, ksi, p, q, r0) / Vect_B_plus.size();

    double E_B_minus = Energy(Vect_B_minus, a0, d_b_minus, A1, A0, ksi, p, q, r0) / Vect_B_minus.size();

    double E_c11_plus = Energy(Vect_c11_plus, a0, d_c11_plus, A1, A0, ksi, p, q, r0) / Vect_c11_plus.size();

    double E_c11_minus = Energy(Vect_c11_minus, a0, d_c11_minus, A1, A0, ksi, p, q, r0) / Vect_c11_minus.size();

    double E_c12_plus = Energy(Vect_c12_plus, a0, d_c12_plus, A1, A0, ksi, p, q, r0) / Vect_c12_plus.size();

    double E_c12_minus = Energy(Vect_c12_minus, a0, d_c12_minus, A1, A0, ksi, p, q, r0) / Vect_c12_minus.size();

    double E_c44_plus = Energy(Vect_c44_plus, a0, d_c44_plus, A1, A0, ksi, p, q, r0) / Vect_c44_plus.size();

    double E_c44_minus = Energy(Vect_c44_minus, a0, d_c44_minus, A1, A0, ksi, p, q, r0) / Vect_c44_minus.size();*/

    double Energy_per_thread[4][9] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    omp_set_dynamic(0);

    #pragma omp parallel num_threads(4)
    {
        int threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();

        int k = Vect_p.size();

        int left_bound = (k*thread_num) / double(threads);
        int right_bound = (k*(thread_num + 1)) / double(threads);

        if(thread_num == threads - 1)
        {
            right_bound = k;
        }

        for(int i=left_bound; i<right_bound; i++)
        {
            Energy_per_thread[thread_num][0] += Energy_parallel(Vect_p, a0, d_edin, A1, A0, ksi, p, q, r0, i) / Vect_p.size();
            Energy_per_thread[thread_num][1] += Energy_parallel(Vect_B_plus, a0, d_b_plus, A1, A0, ksi, p, q, r0, i) / Vect_B_plus.size();
            Energy_per_thread[thread_num][2] += Energy_parallel(Vect_B_minus, a0, d_b_minus, A1, A0, ksi, p, q, r0, i) / Vect_B_minus.size();
            Energy_per_thread[thread_num][3] += Energy_parallel(Vect_c11_plus, a0, d_c11_plus, A1, A0, ksi, p, q, r0, i) / Vect_c11_plus.size();
            Energy_per_thread[thread_num][4] += Energy_parallel(Vect_c11_minus, a0, d_c11_minus, A1, A0, ksi, p, q, r0, i) / Vect_c11_minus.size();
            Energy_per_thread[thread_num][5] += Energy_parallel(Vect_c12_plus, a0, d_c12_plus, A1, A0, ksi, p, q, r0, i) / Vect_c12_plus.size();
            Energy_per_thread[thread_num][6] += Energy_parallel(Vect_c12_minus, a0, d_c12_minus, A1, A0, ksi, p, q, r0, i) / Vect_c12_minus.size();
            Energy_per_thread[thread_num][7] += Energy_parallel(Vect_c44_plus, a0, d_c44_plus, A1, A0, ksi, p, q, r0, i) / Vect_c44_plus.size();
            Energy_per_thread[thread_num][8] += Energy_parallel(Vect_c44_minus, a0, d_c44_minus, A1, A0, ksi, p, q, r0, i) / Vect_c44_minus.size();
        }
    }

    double E_p = 0, E_B_plus = 0, E_B_minus = 0, E_c11_plus = 0, E_c11_minus = 0, E_c12_plus = 0, E_c12_minus = 0, E_c44_plus = 0, E_c44_minus = 0;

    for(int i=0; i<4; i++)
    {
        E_p += Energy_per_thread[i][0];
        E_B_plus += Energy_per_thread[i][1];
        E_B_minus += Energy_per_thread[i][2];
        E_c11_plus += Energy_per_thread[i][3];
        E_c11_minus += Energy_per_thread[i][4];
        E_c12_plus += Energy_per_thread[i][5];
        E_c12_minus += Energy_per_thread[i][6];
        E_c44_plus += Energy_per_thread[i][7];
        E_c44_minus += Energy_per_thread[i][8];
    }

    E = E_p;

    double d2_E_B = (E_B_plus - 2*E_p +E_B_minus)/alpha2;

    double d2_E_c11 = (E_c11_plus - 2*E_p +E_c11_minus)/alpha2;

    double d2_E_c12 =(E_c12_plus - 2*E_p +E_c12_minus)/alpha2;

    double d2_E_c44 =(E_c44_plus - 2*E_p +E_c44_minus)/alpha2;

    B = (d2_E_B*2*const_p)/(9.0*v0);

    c11 = ((d2_E_c11 + d2_E_c12)*const_p)/(2.0*v0);

    c12 = ((d2_E_c11 - d2_E_c12)*const_p)/(2.0*v0);

    c44 = (d2_E_c44*const_p)/(2*v0);
}

double f(double x[], int j, int n)
{
    double E, B, c11, c12, c44;

    Parameters(E, a0_global, B, c11, c12, c44, x[0*(n+1)+j], x[1*(n+1)+j], x[2*(n+1)+j], x[3*(n+1)+j], x[4*(n+1)+j], x[5*(n+1)+j]);

    double sum = 0;

    sum += (a0_global - true_parameters[0])*(a0_global - true_parameters[0]) / (true_parameters[0]*true_parameters[0]);
    sum += (E - true_parameters[1])*(E - true_parameters[1]) / (true_parameters[1]*true_parameters[1]);
    sum += (B - true_parameters[2])*(B - true_parameters[2]) / (true_parameters[2]*true_parameters[2]);
    sum += (c11 - true_parameters[3])*(c11 - true_parameters[3]) / (true_parameters[3]*true_parameters[3]);
    sum += (c12 - true_parameters[4])*(c12 - true_parameters[4]) / (true_parameters[4]*true_parameters[4]);
    sum += (c44 - true_parameters[5])*(c44 - true_parameters[5]) / (true_parameters[5]*true_parameters[5]);

    sum = sqrt(sum/6);

	return sum;
}

double f1(double x[])
{
    double E, B, c11, c12, c44;

    Parameters(E, a0_global, B, c11, c12, c44, x[0], x[1], x[2], x[3], x[4], x[5]);

    double sum = 0;

    sum += (a0_global - true_parameters[0])*(a0_global - true_parameters[0]) / (true_parameters[0]*true_parameters[0]);
    sum += (E - true_parameters[1])*(E - true_parameters[1]) / (true_parameters[1]*true_parameters[1]);
    sum += (B - true_parameters[2])*(B - true_parameters[2]) / (true_parameters[2]*true_parameters[2]);
    sum += (c11 - true_parameters[3])*(c11 - true_parameters[3]) / (true_parameters[3]*true_parameters[3]);
    sum += (c12 - true_parameters[4])*(c12 - true_parameters[4]) / (true_parameters[4]*true_parameters[4]);
    sum += (c44 - true_parameters[5])*(c44 - true_parameters[5]) / (true_parameters[5]*true_parameters[5]);

    sum = sqrt(sum/6);

	return sum;
}

double f_2(double x[], int j, int n)
{
    double sum = 0;

    double E_sol = EnergySol(a0_global, d_edin, x[0*(n+1)+j], x[1*(n+1)+j], x[2*(n+1)+j], x[3*(n+1)+j], x[4*(n+1)+j], x[5*(n+1)+j], Array, E_cohesion);

    sum += (E_sol - true_parameter)*(E_sol - true_parameter) / (true_parameter*true_parameter);

    sum = sqrt(sum);

	return sum;
}

double f1_2(double x[])
{
    double sum = 0;

    double E_sol = EnergySol(a0_global, d_edin, x[0], x[1], x[2], x[3], x[4], x[5], Array, E_cohesion);

    sum += (E_sol - true_parameter)*(E_sol - true_parameter) / (true_parameter*true_parameter);

    sum = sqrt(sum);

	return sum;
}

double f_3(double x[], int j, int n)
{
    double sum = 0;

    double E_In = EnergyIn(a0_global, d_edin, x[0*(n+1)+j], x[1*(n+1)+j], x[2*(n+1)+j], x[3*(n+1)+j], x[4*(n+1)+j], x[5*(n+1)+j], Array2);
    double E_On = EnergyOn(a0_global, d_edin, x[0*(n+1)+j], x[1*(n+1)+j], x[2*(n+1)+j], x[3*(n+1)+j], x[4*(n+1)+j], x[5*(n+1)+j], Array2);

    sum += (E_In - true_params[0])*(E_In - true_params[0]) / (true_params[0]*true_params[0]);
    sum += (E_On - true_params[1])*(E_On - true_params[1]) / (true_params[1]*true_params[1]);

    sum = sqrt(sum/2);

	return sum;
}

double f1_3(double x[])
{
    double sum = 0;

    double E_In = EnergyIn(a0_global, d_edin, x[0], x[1], x[2], x[3], x[4], x[5], Array2);
    double E_On = EnergyOn(a0_global, d_edin, x[0], x[1], x[2], x[3], x[4], x[5], Array2);

    sum += (E_In - true_params[0])*(E_In - true_params[0]) / (true_params[0]*true_params[0]);
    sum += (E_On - true_params[1])*(E_On - true_params[1]) / (true_params[1]*true_params[1]);

    sum = sqrt(sum/2);

	return sum;
}

double* makeSimplex(double x[], int n, double L)
{
    double r1, r2;

    double* S = new double[n*(n+1)];

    r1 = L * ((n - 1.0 + sqrt(n + 1.0)) / (n * sqrt(2.0)));
    r2 = L * ((sqrt(n + 1.0) - 1.0) / (n * sqrt(2.0)));

    for(int i=0; i<n; i++)
    {
        S[i*(n+1)+0] = x[i];
    }

    for(int j=1; j<n+1; j++)
    {
        for(int i=0; i<n; i++)
        {
            if(i == j-1)
            {
                S[i*(n+1)+j] = x[i] + r1;
            }
            else
            {
                S[i*(n+1)+j] = x[i] + r2;
            }
        }
    }

    return S;
}

double* center_of_gravity(double S[], int k, int n)
{
    double* center_of_gravity = new double[n];

    for(int i=0; i<n; i++)
    {
        center_of_gravity[i] = 0.0;
    }

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n+1; j++)
        {
            if(j != k)
            {
                center_of_gravity[i] += S[i*(n+1)+j];
            }
        }

        center_of_gravity[i] = center_of_gravity[i]/n;
    }

    return center_of_gravity;
}

double reflection(double xc[], double xh[], int k, double alpha, int n, int i)
{
    double xr;

    xr = (1.0 + alpha) * xc[i] - alpha * xh[i];

    return xr;
}

double stretch(double xc[], double xr[], int k, double gamma, int n, int i)
{
    double xe;

    xe = (1.0 - gamma) * xc[i] + gamma * xr[i];

    return xe;
}

double compress(double xc[], double xh[], int k, double beta, int n, int i)
{
    double xs;

    xs = (1.0 - beta) * xc[i] + beta * xh[i];

    return xs;
}

bool notStopYet(double (*op1)(double*, int, int), double S[], double eps, int n)
{
    bool notStopYet = false;

    double F[n+1];

    for(int j=0; j<n+1; j++)
    {
        F[j] = op1(S, j, n);
    }

    double min_f = F[0];
    ind_min_global = 0;

    for(int j=0; j<n+1; j++)
    {
        if(F[j] < min_f)
        {
            min_f = F[j];
            ind_min_global = j;
        }
    }

    if(min_f > eps)
        notStopYet = true;

    cout << min_f << endl;

    return notStopYet;
}

double* nelMead(double (*op1)(double*, int, int), double (*op2)(double*), double x[], int n, double L, double eps, double alpha, double beta, double gamma, double sigma)
{
    double *S = makeSimplex(x, n, L);

    double *xc = new double[n];

    double F[n+1], Fr, Fe, Fs;

    double max_f, min_f, max2_f;

    double xh[n], xl[n], xr[n], xe[n], xs[n];

    int ind_max = 0, ind_min;

    bool flag;

    for(int j=0; j<n+1; j++)
    {
        F[j] = op1(S, j, n);
    }

    while(notStopYet(op1, S, eps, n))
    {
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n+1; j++)
            {
                if(S[i*(n+1)+j] < Params_bounds[i][0])
                {
                    S[i*(n+1)+j] = Params_bounds[i][0]*1.01;
                }
                else if(S[i*(n+1)+j] > Params_bounds[i][1])
                {
                    S[i*(n+1)+j] = Params_bounds[i][1]*0.99;
                }
            }
        }

        for(int j=0; j<n+1; j++)
        {
            F[j] = op1(S, j, n);
        }

        max_f = F[0];
        max2_f = F[0];
        min_f = F[0];
        ind_max = 0;
        ind_min = 0;

        for(int j=0; j<n+1; j++)
        {
            if(F[j] <= min_f)
            {
                min_f = F[j];
                ind_min = j;
            }

            if(F[j] >= max_f)
            {
                max2_f = max_f;
                max_f = F[j];
                ind_max = j;
            }
            else if(F[j] >= max2_f)
            {
                max2_f = F[j];
            }
        }

        for(int i=0; i<n; i++)
        {
            xh[i] = S[i*(n+1)+ind_max];
            xl[i] = S[i*(n+1)+ind_min];
        }

        xc = center_of_gravity(S, ind_max, n);

        for(int i=0; i<n; i++)
        {
            xr[i] = reflection(xc, xh, ind_max, alpha, n, i);
        }

        Fr = op2(xr);

        if(Fr < min_f)
        {
            for(int i=0; i<n; i++)
            {
                xe[i] = stretch(xc, xr, ind_max, gamma, n, i);
            }

            Fe = op2(xe);

            if(Fe < Fr)
            {
                for(int i=0; i<n; i++)
                {
                    S[i*(n+1)+ind_max] = xe[i];
                }
                F[ind_max] = Fe;
            }
            else
            {
                for(int i=0; i<n; i++)
                {
                    S[i*(n+1)+ind_max] = xr[i];
                }
                F[ind_max] = Fr;
            }
        }

        if((Fr >= min_f) && (Fr < max2_f))
        {
            for(int i=0; i<n; i++)
            {
                S[i*(n+1)+ind_max] = xr[i];
            }
            F[ind_max] = Fr;
        }

        flag = false;

        if((Fr >= max2_f) && (Fr < max_f))
        {
            flag = true;

            for(int i=0; i<n; i++)
            {
                S[i*(n+1)+ind_max] = xr[i];
            }
            F[ind_max] = Fr;

            for(int i=0; i<n; i++)
            {
                xh[i] = xr[i];
            }
        }

        if(Fr >= max_f)
        {
            flag = true;
        }

        if(flag)
        {
            for(int i=0; i<n; i++)
            {
                xs[i] = compress(xc, xh, ind_max, beta, n, i);
            }
            Fs = op2(xs);

            if(Fs < max_f)
            {
                for(int i=0; i<n; i++)
                {
                    S[i*(n+1)+ind_max] = xs[i];
                }
                F[ind_max] = Fs;
            }
            else
            {
                for(int j=0; j<n+1; j++)
                {
                    for(int i=0; i<n; i++)
                    {
                        S[i*(n+1)+j] = xl[i] + sigma*(S[i*(n+1)+j] - xl[i]);
                    }
                }

                for(int j=0; j<n+1; j++)
                {
                    F[j] = op1(S, j, n);
                }
            }
        }
    }

    for(int j=0; j<n+1; j++)
    {
        F[j] = op1(S, j, n);
    }

    double min_f_global = F[0];
    ind_min_global = 0;

    for(int j=0; j<n+1; j++)
    {
        if(F[j] < min_f_global)
        {
            min_f_global = F[j];
            ind_min_global = j;
        }
    }

    double* Params = new double[n];

    for(int i=0; i<n; i++)
    {
        Params[i] = S[i*(n+1)+ind_min_global];
    }

    return Params;
}

double GetRandomNumber(double min, double max, int precision)
{
  double value;

  value = rand() % (int)pow(10, precision);

  value = min + (value / pow(10, precision)) * (max - min);

  return value;
}

int main()
{
    srand(4259707270);

    int n = 6;

    cout << endl;
    cout << "   A1: " << Params_bounds[0][0] << " - " << Params_bounds[0][1] << endl;
    cout << "   A0: " << Params_bounds[1][0] << " - " << Params_bounds[1][1] << endl;
    cout << "   ksi: " << Params_bounds[2][0] << " - " << Params_bounds[2][1] << endl;
    cout << "   p: " << Params_bounds[3][0] << " - " << Params_bounds[3][1] << endl;
    cout << "   q: " << Params_bounds[4][0] << " - " << Params_bounds[4][1] << endl;
    cout << "   r0: " << Params_bounds[5][0] << " - " << Params_bounds[5][1] << endl;

    ofstream fout1;
    ofstream fout2;
    ofstream fout3;
    ofstream fout4;

	fout1.open("x.txt");
	fout2.open("y1.txt");
	fout3.open("y2.txt");
	fout4.open("y3.txt");

	double* Params1 = new double[n];
	double* Params2 = new double[n];
	double* Params3 = new double[n];

    a0_global = rand()%11+2;

    h = 1.0;

	double x[n];

	unsigned int start_time =  clock();

	x[0] = 0.0280413;
	x[1] = 0.110633;
	x[2] = 1.19174;
	x[3] = 11.1364;
	x[4] = 3.03011;
	x[5] = 2.86881;

	while(h > 1.0e-6)
    {
        double a_left, a_right;
        double Energy_left, Energy_right, Energy_cent;

        a_left = a0_global - h;
        a_right = a0_global + h;

        vector <Atom> Vect_left = GCK("Ag", a_left);
        vector <Atom> Vect_right = GCK("Ag", a_right);
        vector <Atom> Vect_cent = GCK("Ag", a0_global);

        Energy_left = Energy(Vect_left, a_left, d_edin, x[0], x[1], x[2], x[3], x[4], x[5]);
        Energy_right = Energy(Vect_right, a_right, d_edin, x[0], x[1], x[2], x[3], x[4], x[5]);
        Energy_cent = Energy(Vect_cent, a0_global, d_edin, x[0], x[1], x[2], x[3], x[4], x[5]);

        double Energy_min = min(min(Energy_left, Energy_right), Energy_cent);

        if(Energy_min == Energy_cent)
        {
            a0_global = a0_global;
            h /= 10;
        }
        else if(Energy_min == Energy_left)
        {
            a0_global = a_left;
        }
        else if(Energy_min == Energy_right)
        {
            a0_global = a_right;
        }
    }

	x[0] = GetRandomNumber(0.0, 0.1, 6);
	x[1] = GetRandomNumber(0.0685, 0.1370, 6);
	x[2] = GetRandomNumber(0.7853, 1.570666, 6);
	x[3] = GetRandomNumber(7.2853, 14.5706, 6);
	x[4] = GetRandomNumber(2.0927, 4.1853, 6);
	x[5] = GetRandomNumber(1.9257, 3.8514, 6);

    Params1 = nelMead(f, f1, x, n, 1.0, 0.003, 1.0, 0.5, 2.0, 0.5);

    cout << endl;
    cout << "   B-B:" << endl;
    cout << "   A1 = " << Params1[0] << endl;
    cout << "   A0 = " << Params1[1] << endl;
    cout << "   ksi = " << Params1[2] << endl;
    cout << "   p = " << Params1[3] << endl;
    cout << "   q = " << Params1[4] << endl;
    cout << "   r0 = " << Params1[5] << endl;

    double E, a0, B, c11, c12, c44;

    a0 = a0_global;

    Parameters(E, a0_global, B, c11, c12, c44, Params1[0], Params1[1], Params1[2], Params1[3], Params1[4], Params1[5]);

    cout << endl;
    cout << "   a0 = " << a0 << endl  << "   E_coh = " << E << endl << "   B = " << B << endl << "   c11 = " << c11 << endl << "   c12 = " << c12 << endl << "   c44 = " << c44 << endl;
    cout << endl;

    for(int i=0; i<n; i++)
    {
        Array[i] = Params1[i];
    }
    E_cohesion = E;

    x[0] = GetRandomNumber(0.0, 0.1, 6);
	x[1] = GetRandomNumber(0.0685, 0.1370, 6);
	x[2] = GetRandomNumber(0.7853, 1.570666, 6);
	x[3] = GetRandomNumber(7.2853, 14.5706, 6);
	x[4] = GetRandomNumber(2.0927, 4.1853, 6);
	x[5] = GetRandomNumber(1.9257, 3.8514, 6);

    Params2 = nelMead(f_2, f1_2, x, n, 1.0, 0.000001, 1.0, 0.5, 2.0, 0.5);

    cout << "   A-B:" << endl;
    cout << "   A1 = " << Params2[0] << endl;
    cout << "   A0 = " << Params2[1] << endl;
    cout << "   ksi = " << Params2[2] << endl;
    cout << "   p = " << Params2[3] << endl;
    cout << "   q = " << Params2[4] << endl;
    cout << "   r0 = " << Params2[5] << endl;

    cout << endl;
    cout << "   E_sol = " << EnergySol(a0_global, d_edin, Params2[0], Params2[1], Params2[2], Params2[3], Params2[4], Params2[5], Array, E_cohesion) << endl;
    cout << endl;

    for(int i=0; i<n; i++)
    {
        Array2[i] = Params1[i];
    }

    for(int i=n; i<2*n; i++)
    {
        Array2[i] = Params2[i-n];
    }

    x[0] = GetRandomNumber(0.0, 0.1, 6);
	x[1] = GetRandomNumber(0.0685, 0.1370, 6);
	x[2] = GetRandomNumber(0.7853, 1.570666, 6);
	x[3] = GetRandomNumber(7.2853, 14.5706, 6);
	x[4] = GetRandomNumber(2.0927, 4.1853, 6);
	x[5] = GetRandomNumber(1.9257, 3.8514, 6);

    Params3 = nelMead(f_3, f1_3, x, n, 1.0, 0.000001, 1.0, 0.5, 2.0, 0.5);

    cout << "   A-A:" << endl;
    cout << "   A1 = " << Params3[0] << endl;
    cout << "   A0 = " << Params3[1] << endl;
    cout << "   ksi = " << Params3[2] << endl;
    cout << "   p = " << Params3[3] << endl;
    cout << "   q = " << Params3[4] << endl;
    cout << "   r0 = " << Params3[5] << endl;

    cout << endl;
    cout << "   E_in = " << EnergyIn(a0_global, d_edin, Params3[0], Params3[1], Params3[2], Params3[3], Params3[4], Params3[5], Array2) << endl;
    cout << "   E_on = " << EnergyOn(a0_global, d_edin, Params3[0], Params3[1], Params3[2], Params3[3], Params3[4], Params3[5], Array2) << endl;
    cout << endl;

    unsigned int end_time = clock();

    unsigned int search_time = end_time - start_time;
    search_time = search_time / 1000.0;

    cout << "   Execution time: " << search_time << " seconds" << endl;
    cout << endl;

    double rij = 0.25;

    while(rij <= 6.5)
    {
        fout1 << rij << endl;

        fout2 << EnergyFinal(rij, Params1) << endl;
        fout3 << EnergyFinal(rij, Params2) << endl;
        fout4 << EnergyFinal(rij, Params3) << endl;

        rij += 0.05;
    }

    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
}
